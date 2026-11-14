#include "ramcore/RecordBlocks.h"
#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <htscodecs/fqzcomp_qual.h>
#include <htscodecs/rANS_static4x16.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <future>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

// After <cstdint>: the header uses uint8_t without including it.
#include <htscodecs/tokenise_name3.h>

namespace {

// CRAM 3.1's variant: reverse-strand qualities are modelled in sequencing order.
constexpr int kFqzVersion = 3;
// fqzcomp's parameter set for MiSeq-like data; the smallest on HG00154.
constexpr int kFqzStrategy = 2;
// The name tokeniser's level, of 1 to 9. On NA12878 chr20, level 9 with
// arithmetic coding saves 2 MB more of 61 MB, but encodes 4.7 times slower.
constexpr int kTok3Level = 3;
// rANS Nx16 order 1 with RLE and PACK, as CRAM 3.1 compresses its quality
// series. On NA12878 chr20: 72.4 MB, against 71.6 MB with fqzcomp, which
// decodes 25 times slower, and 114.6 MB for the qualities as a zstd column.
constexpr int kRansOrder = 1 | RANS_ORDER_RLE | RANS_ORDER_PACK;

// The first byte of a quality block says how the rest is compressed. A rANS
// block holds the size of its lengths part in four bytes, little endian, the
// lengths as varints, then the qualities, each part compressed with kRansOrder.
// An fqzcomp block is fqzcomp's own stream, which holds the lengths.
constexpr char kQualRans = 1;
constexpr char kQualFqzcomp = 2;

bool IsSamQuality(const std::string &qual)
{
   return std::all_of(qual.begin(), qual.end(), [](char c) { return c >= '!' && c <= '~'; });
}

// fqzcomp and the name tokeniser take char, rANS unsigned char.
unsigned char *AsBytes(char *p)
{
   return static_cast<unsigned char *>(static_cast<void *>(p));
}

std::string Rans(std::string &in)
{
   unsigned int size = 0;
   unsigned char *out = rans_compress_4x16(AsBytes(in.data()), static_cast<unsigned int>(in.size()), &size, kRansOrder);
   if (!out)
      throw std::runtime_error("rANS could not compress a quality block");
   std::string bytes(size, '\0');
   std::memcpy(bytes.data(), out, size);
   free(out); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs
   return bytes;
}

std::string Unrans(char *in, std::size_t size)
{
   unsigned int n = 0;
   unsigned char *out = rans_uncompress_4x16(AsBytes(in), static_cast<unsigned int>(size), &n);
   if (!out)
      throw std::runtime_error("rANS could not decode a quality block");
   std::string bytes(n, '\0');
   std::memcpy(bytes.data(), out, n);
   free(out); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs
   return bytes;
}

// The block of \p lengths and \p quals, as kQualRans describes it.
std::string RansBlock(const std::vector<uint32_t> &lengths, std::string &quals)
{
   std::string varints;
   for (uint32_t n : lengths) {
      for (; n >= 0x80; n >>= 7)
         varints.push_back(static_cast<char>((n & 0x7f) | 0x80));
      varints.push_back(static_cast<char>(n));
   }
   const std::string packedLengths = Rans(varints);
   std::string block(1, kQualRans);
   for (int shift = 0; shift < 32; shift += 8)
      block.push_back(static_cast<char>((packedLengths.size() >> shift) & 0xff));
   block += packedLengths;
   block += Rans(quals);
   return block;
}

// Fills \p lengths and \p quals from the block after its codec byte, for \p records qualities.
void DecodeRansBlock(std::vector<char> &bytes, std::size_t records, std::vector<int> &lengths, std::string &quals)
{
   constexpr std::size_t kHeader = 5;
   if (bytes.size() < kHeader)
      throw std::runtime_error("a quality block is damaged");
   std::size_t lengthsSize = 0;
   for (std::size_t i = 0; i < 4; i++)
      lengthsSize |= static_cast<std::size_t>(static_cast<unsigned char>(bytes.at(1 + i))) << (8 * i);
   if (lengthsSize >= bytes.size() - kHeader)
      throw std::runtime_error("a quality block is damaged");

   const std::string varints = Unrans(&bytes.at(kHeader), lengthsSize);
   lengths.clear();
   std::size_t total = 0;
   uint32_t n = 0;
   int shift = 0;
   for (const char c : varints) {
      const auto byte = static_cast<unsigned char>(c);
      if (shift > 28)
         throw std::runtime_error("a quality block is damaged");
      n |= static_cast<uint32_t>(byte & 0x7f) << shift;
      shift += 7;
      if (byte < 0x80) {
         if (n > static_cast<uint32_t>(std::numeric_limits<int>::max()))
            throw std::runtime_error("a quality block is damaged");
         lengths.push_back(static_cast<int>(n));
         total += n;
         n = 0;
         shift = 0;
      }
   }
   quals = Unrans(&bytes.at(kHeader + lengthsSize), bytes.size() - kHeader - lengthsSize);
   if (shift != 0 || lengths.size() != records || quals.size() != total)
      throw std::runtime_error("a quality block is damaged");
}

// Ends that do not rise, or that run past the last row, would send a row to a
// block that does not hold it.
void CheckBlockEnds(const std::vector<uint64_t> &ends, ROOT::NTupleSize_t rows)
{
   for (std::size_t i = 0; i < ends.size(); i++) {
      if ((i > 0 && ends[i] <= ends[i - 1]) || ends[i] >= rows)
         throw std::runtime_error("the block ends in the file are damaged");
   }
}

// What a CIGAR operation consumes, as htslib's bam_cigar_type() says it.
constexpr int kRead = 1;
constexpr int kRef = 2;
constexpr int kReadAndRef = kRead | kRef;

// kRead, kRef, both or neither; -1 for an operation SAM does not have.
int CigarType(uint32_t op)
{
   switch (op & 0xF) {
   case RAM_CIGAR_M:
   case RAM_CIGAR_EQUAL:
   case RAM_CIGAR_X: return kReadAndRef;
   case RAM_CIGAR_I:
   case RAM_CIGAR_S: return kRead;
   case RAM_CIGAR_D:
   case RAM_CIGAR_N: return kRef;
   case RAM_CIGAR_H:
   case RAM_CIGAR_P: return 0;
   default: return -1;
   }
}

constexpr std::string_view kAcgt = "ACGT";
// SEQ's alphabet in htslib's 4-bit order; RAMNTupleUtils::NormalizeSequence() gives only these.
constexpr std::string_view kBases = "=ACMGRSVTWYHKDBN";
// A seq_diff code below kSubOther picks one of the three bases other than the
// consensus's A, C, G or T, as CRAM codes substitutions; from kSubOther up it
// is kBases[code - kSubOther].
constexpr int kSubOther = 3;
// A seq_diff holds its code in the low bits and its step from the previous
// difference above them.
constexpr unsigned kDiffCodeBits = 5;
constexpr std::size_t kMaxDiffStep = (std::size_t{1} << (32 - kDiffCodeBits)) - 1;

// A, C, G and T as 0 to 3, every other byte as -1: a table, since the writer
// looks up every aligned base.
constexpr std::array<std::int8_t, 256> kAcgtIndex = [] {
   std::array<std::int8_t, 256> t{};
   for (auto &v : t)
      v = -1;
   for (std::size_t i = 0; i < kAcgt.size(); i++)
      t[static_cast<unsigned char>(kAcgt[i])] = static_cast<std::int8_t>(i);
   return t;
}();

// A count's slot for each byte: A, C, G and T, then one for every other byte,
// which no consensus base comes from. It spares the counting loop a branch.
constexpr std::size_t kOtherSlot = 4;
constexpr std::array<std::uint8_t, 256> kCountSlot = [] {
   std::array<std::uint8_t, 256> t{};
   for (auto &v : t)
      v = kOtherSlot;
   for (std::size_t i = 0; i < kAcgt.size(); i++)
      t[static_cast<unsigned char>(kAcgt[i])] = static_cast<std::uint8_t>(i);
   return t;
}();

// The base each seq_diff code gives in place of each byte of the consensus, or
// '\0' where it gives none: a table, since a long read has thousands of
// differences to apply.
constexpr std::array<std::array<char, 256>, std::size_t{1} << kDiffCodeBits> kSubstitute = [] {
   std::array<std::array<char, 256>, std::size_t{1} << kDiffCodeBits> t{};
   for (std::size_t code = 0; code < static_cast<std::size_t>(kSubOther); code++) {
      for (std::size_t cons = 0; cons < kAcgt.size(); cons++)
         t[code][static_cast<unsigned char>(kAcgt[cons])] = kAcgt[(cons + code + 1) % 4];
   }
   for (std::size_t i = 0; i < kBases.size(); i++) {
      for (auto &base : t[static_cast<std::size_t>(kSubOther) + i])
         base = kBases[i];
   }
   return t;
}();

// The operations CigarType() calls kReadAndRef, a bit each: a test of one is
// cheaper than its switch, in a loop over thousands of them a long read.
constexpr uint32_t kReadAndRefOps = (1U << RAM_CIGAR_M) | (1U << RAM_CIGAR_EQUAL) | (1U << RAM_CIGAR_X);

int Acgt(char base)
{
   return kAcgtIndex[static_cast<unsigned char>(base)];
}

// Whether the SEQ of \p rec can be coded against a consensus: a placed, mapped
// read whose CIGAR accounts for every base. If so, sets \p end to its reference end.
bool Codable(const RAMNTupleRecord &rec, std::int64_t &end)
{
   if (rec.refid < 0 || rec.pos < 0 || (rec.flag & 0x4) || rec.cigar.empty() || rec.seq.empty() || rec.seq == "*")
      return false;
   std::uint64_t bases = 0;
   std::int64_t span = 0;
   for (const uint32_t op : rec.cigar) {
      const int type = CigarType(op);
      if (type < 0)
         return false;
      if (type & kRead)
         bases += op >> 4;
      if (type & kRef)
         span += op >> 4;
   }
   if (bases != rec.seq.size())
      return false;
   end = rec.pos + span;
   return true;
}

// The seq_diff code of \p base where the consensus has \p cons; -1 for a base outside kBases.
int SubCode(char base, char cons)
{
   const int b = Acgt(base);
   const int c = Acgt(cons);
   if (b >= 0 && c >= 0)
      return (b - c + 4) % 4 - 1;
   const auto i = kBases.find(base);
   return i == std::string_view::npos ? -1 : kSubOther + static_cast<int>(i);
}

// Fills the stored SEQ of \p rec, which Codable() accepted, against \p consensus
// starting at \p lo. False, leaving the fields half filled, for a base outside kBases.
bool CodeAgainst(RAMNTupleRecord &rec, std::string_view consensus, std::int64_t lo)
{
   std::size_t q = 0;
   auto at = static_cast<std::size_t>(rec.pos - lo);
   std::size_t last = 0;
   for (const uint32_t op : rec.cigar) {
      const int type = CigarType(op);
      const std::size_t n = op >> 4;
      if (type == kReadAndRef) {
         // Most bases match: skip to each difference, 16 bytes at a time.
         constexpr std::size_t kStride = 16;
         const char *read = rec.seq.data() + q;
         const char *cons = consensus.data() + at;
         for (std::size_t k = 0; k < n; k++) {
            while (n - k >= kStride && std::memcmp(read + k, cons + k, kStride) == 0)
               k += kStride;
            while (k < n && read[k] == cons[k])
               k++;
            if (k == n)
               break;
            const int code = SubCode(read[k], cons[k]);
            if (code < 0)
               return false;
            if (q + k - last > kMaxDiffStep)
               return false;
            rec.seq_diff.push_back(static_cast<uint32_t>((q + k - last) << kDiffCodeBits) |
                                   static_cast<uint32_t>(code));
            last = q + k;
         }
      } else if (type & kRead) {
         rec.seq_lit.append(rec.seq, q, n);
      }
      if (type & kRead)
         q += n;
      if (type & kRef)
         at += n;
   }
   return true;
}

} // namespace

BlockWriter::BlockWriter(std::unique_ptr<ROOT::REntry> entry, FillFn fill, std::size_t block_records,
                         std::size_t block_bases, bool background)
   : m_entry(std::move(entry)),
     m_recordToken(m_entry->GetToken("record")),
     m_qualBlob(m_entry->GetPtr<std::vector<std::uint8_t>>(RAMNTupleRecord::kQualBlockField).get()),
     m_nameBlob(m_entry->GetPtr<std::vector<std::uint8_t>>(RAMNTupleRecord::kNameBlockField).get()),
     m_consensus(m_entry->GetPtr<std::string>(RAMNTupleRecord::kConsensusField).get()),
     m_consensusPos(m_entry->GetPtr<std::int32_t>(RAMNTupleRecord::kConsensusPosField).get()),
     m_fill(std::move(fill)),
     m_blockRecords(std::max<std::size_t>(block_records, 1)),
     m_blockBases(block_bases),
     m_background(background)
{
   // Fill() binds records untyped, which skips ROOT's type check per row (27
   // string compares a record); check the type once here instead.
   m_entry->GetPtr<RAMNTupleRecord>(m_recordToken);
   m_open.records.reserve(std::min(m_blockRecords, kBlockRecords));
}

BlockWriter::~BlockWriter()
{
   // The encoding thread uses m_closed; wait for it, and drop what it threw:
   // a destructor runs with a block still encoding only when the conversion
   // has already failed.
   if (m_encoding.valid())
      m_encoding.wait();
}

RAMNTupleRecord &BlockWriter::Record()
{
   if (m_open.records.size() == m_open.rows)
      m_open.records.emplace_back();
   return m_open.records[m_open.rows];
}

void BlockWriter::Add()
{
   RAMNTupleRecord &rec = m_open.records.at(m_open.rows);
   m_open.names += rec.qname;
   m_open.names += '\0';

   // kFqzcomp comes in from the quality policy and picks the block's codec;
   // kQualInBlock goes out meaning the quality is in the block. fqzcomp takes
   // no empty strings, so "*" stays in the record. So do the lossy policies:
   // their bytes only read back through DecodeQuality, and binned qualities
   // from Q30 up look like SAM quality text.
   if (rec.TestBit(RAMNTupleRecord::kFqzcomp))
      m_open.fqzcomp = true;
   rec.compression_flags &= ~static_cast<uint32_t>(RAMNTupleRecord::kFqzcomp);
   const bool lossless = !rec.TestBit(RAMNTupleRecord::kIlluminaBinning | RAMNTupleRecord::kDrop);
   if (lossless && rec.qual != "*" && !rec.qual.empty() && IsSamQuality(rec.qual)) {
      for (const char c : rec.qual)
         m_open.quals.push_back(static_cast<char>(c - '!'));
      m_open.lengths.push_back(static_cast<uint32_t>(rec.qual.size()));
      m_open.flags.push_back(((rec.flag & 0x10) ? FQZ_FREVERSE : 0) | ((rec.flag & 0x80) ? FQZ_FREAD2 : 0));
      rec.qual.clear();
      rec.SetBit(RAMNTupleRecord::kQualInBlock);
   } else {
      rec.compression_flags &= ~static_cast<uint32_t>(RAMNTupleRecord::kQualInBlock);
   }
   m_open.bases += rec.seq.size();
   m_open.rows++;
   if (m_open.rows == m_blockRecords || m_open.quals.size() >= m_blockBases || m_open.bases >= m_blockBases)
      EndBlock();
}

void BlockWriter::Finish()
{
   if (m_open.rows > 0)
      EndBlock();
   if (m_encoding.valid()) {
      m_encoding.get();
      Fill(m_closed);
   }
}

void BlockWriter::EndBlock()
{
   if (m_encoding.valid()) {
      m_encoding.get();
      Fill(m_closed);
   }
   // Swapping keeps both blocks' buffers for reuse.
   std::swap(m_open, m_closed);
   if (m_background) {
      m_encoding = std::async(std::launch::async, [this] { Encode(m_closed); });
   } else {
      Encode(m_closed);
      Fill(m_closed);
   }
}

void BlockWriter::Encode(Block &b)
{
   CodeSequences(b);
   if (!b.lengths.empty() && b.fqzcomp) {
      fqz_slice slice{};
      slice.num_records = static_cast<int>(b.lengths.size());
      slice.len = b.lengths.data();
      slice.flags = b.flags.data();
      std::size_t size = 0;
      char *packed =
         fqz_compress(kFqzVersion, &slice, b.quals.data(), b.quals.size(), &size, kFqzStrategy, /*gp=*/nullptr);
      if (!packed)
         throw std::runtime_error("fqzcomp could not compress a quality block");
      const std::string_view bytes(packed, size);
      b.qualBlob.assign(1, kQualFqzcomp);
      b.qualBlob.insert(b.qualBlob.end(), bytes.begin(), bytes.end());
      free(packed); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs
   } else if (!b.lengths.empty()) {
      const std::string bytes = RansBlock(b.lengths, b.quals);
      b.qualBlob.assign(bytes.begin(), bytes.end());
   }
   int size = 0;
   int lastStart = 0;
   uint8_t *names = tok3_encode_names(b.names.data(), static_cast<int>(b.names.size()), kTok3Level,
                                      /*use_arith=*/0, &size, &lastStart);
   if (!names)
      throw std::runtime_error("the name tokeniser could not compress a block of names");
   b.nameBlob.resize(static_cast<std::size_t>(size));
   std::memcpy(b.nameBlob.data(), names, static_cast<std::size_t>(size));
   free(names); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs
}

void BlockWriter::Fill(Block &b)
{
   // RNTuple copies a record's values when it fills a row, so binding the
   // block's records in turn saves copying them into the entry.
   for (std::size_t i = 0; i + 1 < b.rows; i++) {
      m_entry->BindRawPtr<void>(m_recordToken, &b.records[i]);
      m_fill(*m_entry);
   }
   m_qualBlob->swap(b.qualBlob);
   m_nameBlob->swap(b.nameBlob);
   m_consensus->swap(b.consensus);
   *m_consensusPos = b.consensusPos;
   m_entry->BindRawPtr<void>(m_recordToken, &b.records[b.rows - 1]);
   m_fill(*m_entry);
   m_qualBlob->clear();
   m_nameBlob->clear();
   m_consensus->clear();
   *m_consensusPos = 0;

   m_rowsSinceTake += b.rows;
   m_blockEnds.push_back(m_rowsSinceTake - 1);
   b.rows = 0;
   b.bases = 0;
   b.names.clear();
   b.quals.clear();
   b.lengths.clear();
   b.flags.clear();
   b.fqzcomp = false;
   b.qualBlob.clear();
   b.nameBlob.clear();
}

void BlockWriter::CodeSequences(Block &b)
{
   // Reads are coded against the reference that has the most of the block's
   // aligned bases; reads on the others keep SEQ as text, as do reads past
   // kMaxConsensusSpan. ends[i] is the reference end of a read to code, else -1.
   std::vector<std::int64_t> ends(b.rows, -1);
   std::vector<std::pair<std::int32_t, std::uint64_t>> basesByRef;
   for (std::size_t i = 0; i < b.rows; i++) {
      const RAMNTupleRecord &rec = b.records[i];
      if (!Codable(rec, ends[i]))
         continue;
      auto it = std::find_if(basesByRef.begin(), basesByRef.end(), [&](const auto &r) { return r.first == rec.refid; });
      if (it == basesByRef.end())
         it = basesByRef.insert(basesByRef.end(), {rec.refid, 0});
      it->second += rec.seq.size();
   }
   const auto most = std::max_element(basesByRef.begin(), basesByRef.end(),
                                      [](const auto &x, const auto &y) { return x.second < y.second; });
   const std::int32_t refid = most == basesByRef.end() ? -1 : most->first;

   std::int64_t lo = std::numeric_limits<std::int64_t>::max();
   for (std::size_t i = 0; i < b.rows; i++) {
      if (ends[i] >= 0 && b.records[i].refid != refid)
         ends[i] = -1;
      if (ends[i] >= 0)
         lo = std::min<std::int64_t>(lo, b.records[i].pos);
   }
   std::int64_t hi = lo;
   for (std::size_t i = 0; i < b.rows; i++) {
      if (ends[i] >= 0 && ends[i] - lo > kMaxConsensusSpan)
         ends[i] = -1;
      if (ends[i] >= 0)
         hi = std::max(hi, ends[i]);
   }

   // The consensus: the commonest of A, C, G and T at each base, N where none.
   std::string &consensus = b.consensus;
   consensus.clear();
   if (refid >= 0) {
      b.counts.assign(static_cast<std::size_t>(hi - lo), {});
      for (std::size_t i = 0; i < b.rows; i++) {
         if (ends[i] < 0)
            continue;
         const RAMNTupleRecord &rec = b.records[i];
         std::size_t q = 0;
         auto at = static_cast<std::size_t>(rec.pos - lo);
         for (const uint32_t op : rec.cigar) {
            const int type = CigarType(op);
            const std::size_t n = op >> 4;
            if (type == kReadAndRef) {
               const char *read = rec.seq.data() + q;
               auto *count = b.counts.data() + at;
               for (std::size_t k = 0; k < n; k++)
                  count[k][kCountSlot[static_cast<unsigned char>(read[k])]]++;
            }
            if (type & kRead)
               q += n;
            if (type & kRef)
               at += n;
         }
      }
      consensus.assign(b.counts.size(), 'N');
      for (std::size_t i = 0; i < b.counts.size(); i++) {
         const auto &c = b.counts[i];
         const auto best = std::max_element(c.begin(), c.begin() + kOtherSlot);
         if (*best > 0)
            consensus[i] = kAcgt[static_cast<std::size_t>(best - c.begin())];
      }
   }

   for (std::size_t i = 0; i < b.rows; i++) {
      RAMNTupleRecord &rec = b.records[i];
      rec.seq_lit.clear();
      rec.seq_diff.clear();
      rec.compression_flags &= ~static_cast<uint32_t>(RAMNTupleRecord::kSeqConsensus);
      if (ends[i] >= 0 && CodeAgainst(rec, consensus, lo)) {
         rec.SetBit(RAMNTupleRecord::kSeqConsensus);
      } else {
         rec.seq_lit = rec.seq;
         rec.seq_diff.clear();
      }
   }
   b.consensusPos = refid < 0 ? 0 : static_cast<std::int32_t>(lo);
}

std::vector<uint64_t> BlockWriter::TakeBlockEnds()
{
   m_rowsSinceTake = 0;
   return std::exchange(m_blockEnds, {});
}

QualityBlockReader::QualityBlockReader(RAMFile &file)
   : m_flagsView(file.GetReader().GetView<uint32_t>("record.compression_flags")),
     m_ends(file.GetMetadata().block_ends),
     m_readAhead(std::min<std::size_t>(kQualityReadAhead, std::max(1U, std::thread::hardware_concurrency())))
{
   ROOT::RNTupleReader &reader = file.GetReader();
   if (reader.GetDescriptor().FindFieldId(RAMNTupleRecord::kQualBlockField) != ROOT::kInvalidDescriptorId)
      m_view.emplace(reader.GetView<std::vector<std::uint8_t>>(RAMNTupleRecord::kQualBlockField));
   CheckBlockEnds(m_ends, reader.GetNEntries());
}

std::string QualityBlockReader::Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row)
{
   if (!rec.TestBit(RAMNTupleRecord::kQualInBlock))
      return rec.GetQUAL();

   const auto it = std::lower_bound(m_ends.begin(), m_ends.end(), static_cast<uint64_t>(row));
   if (it == m_ends.end() || !m_view)
      throw std::runtime_error("record " + std::to_string(row) + " has no quality block");
   const auto block = static_cast<std::size_t>(it - m_ends.begin());
   if (block != m_block)
      Load(block);

   const std::size_t slot = m_current.slots[static_cast<std::size_t>(row - m_firstRow)];
   return m_current.quals.substr(m_current.offsets[slot], static_cast<std::size_t>(m_current.lengths[slot]));
}

QualityBlockReader::Packed QualityBlockReader::Read(std::size_t block)
{
   const uint64_t first = block == 0 ? 0 : m_ends[block - 1] + 1;
   Packed p;
   p.slots.assign(static_cast<std::size_t>(m_ends[block] - first + 1), 0);
   for (std::size_t i = 0; i < p.slots.size(); i++) {
      if (m_flagsView(first + i) & RAMNTupleRecord::kQualInBlock)
         p.slots[i] = p.records++;
   }
   if (p.records > 0) {
      const std::vector<std::uint8_t> &bytes = (*m_view)(m_ends[block]);
      p.bytes.assign(bytes.begin(), bytes.end());
   }
   return p;
}

QualityBlockReader::Block QualityBlockReader::Decode(Packed packed)
{
   Block b;
   b.slots = std::move(packed.slots);
   b.lengths.assign(packed.records, 0);
   if (packed.records > 0) {
      std::vector<char> &bytes = packed.bytes;
      if (bytes.size() < 2)
         throw std::runtime_error("a quality block is damaged");
      if (bytes.front() == kQualRans) {
         DecodeRansBlock(bytes, packed.records, b.lengths, b.quals);
      } else if (bytes.front() == kQualFqzcomp) {
         std::size_t size = 0;
         char *out =
            fqz_decompress(&bytes.at(1), bytes.size() - 1, &size, b.lengths.data(), static_cast<int>(packed.records));
         if (!out)
            throw std::runtime_error("fqzcomp could not decode a quality block");
         b.quals.assign(out, size);
         free(out); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs
      } else {
         throw std::runtime_error("a quality block has an unknown codec");
      }
      for (auto &c : b.quals)
         c = static_cast<char>(c + '!');
   }
   b.offsets.resize(packed.records);
   std::size_t offset = 0;
   for (std::size_t i = 0; i < packed.records; i++) {
      b.offsets[i] = offset;
      offset += static_cast<std::size_t>(b.lengths[i]);
   }
   return b;
}

void QualityBlockReader::Load(std::size_t block)
{
   // Read ahead one block more for each consecutive block, so random lookups start nothing.
   m_run = (m_block != kNoBlock && block == m_block + 1) ? m_run + 1 : 0;

   auto ahead = m_ahead.find(block);
   if (ahead != m_ahead.end()) {
      m_current = ahead->second.get();
      m_ahead.erase(ahead);
   } else {
      m_current = Decode(Read(block));
   }
   m_block = block;
   m_firstRow = block == 0 ? 0 : m_ends[block - 1] + 1;

   const std::size_t depth = std::min(m_run, m_readAhead);
   for (std::size_t next = block + 1; next < m_ends.size() && next <= block + depth; next++) {
      if (m_ahead.count(next) == 0)
         m_ahead.emplace(next, std::async(std::launch::async, Decode, Read(next)));
   }
}

NameBlockReader::NameBlockReader(RAMFile &file)
   : m_view(file.GetReader().GetView<std::vector<std::uint8_t>>(RAMNTupleRecord::kNameBlockField)),
     m_ends(file.GetMetadata().block_ends)
{
   CheckBlockEnds(m_ends, file.GetReader().GetNEntries());
}

std::string_view NameBlockReader::Get(ROOT::NTupleSize_t row)
{
   if (m_block == kNoBlock || row < m_firstRow || row > m_ends[m_block]) {
      const auto it = std::lower_bound(m_ends.begin(), m_ends.end(), static_cast<uint64_t>(row));
      if (it == m_ends.end())
         throw std::runtime_error("record " + std::to_string(row) + " has no block of names");
      Load(static_cast<std::size_t>(it - m_ends.begin()));
   }
   const auto i = static_cast<std::size_t>(row - m_firstRow);
   const std::vector<std::size_t> &starts = m_current.starts;
   return std::string_view(m_current.names).substr(starts[i], starts[i + 1] - starts[i] - 1);
}

std::vector<std::uint8_t> NameBlockReader::Read(std::size_t block)
{
   const std::vector<std::uint8_t> &bytes = m_view(m_ends[block]);
   return {bytes.begin(), bytes.end()};
}

NameBlockReader::Block NameBlockReader::Decode(std::vector<std::uint8_t> packed, std::size_t block, std::size_t rows)
{
   uint32_t size = 0;
   uint8_t *names = tok3_decode_names(packed.data(), static_cast<uint32_t>(packed.size()), &size);
   if (!names)
      throw std::runtime_error("the name tokeniser could not decode block " + std::to_string(block));
   Block b;
   b.names.resize(size);
   std::memcpy(b.names.data(), names, size);
   free(names); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs

   b.starts.reserve(rows + 1);
   b.starts.push_back(0);
   for (std::size_t end = b.names.find('\0'); end != std::string::npos; end = b.names.find('\0', end + 1))
      b.starts.push_back(end + 1);
   if (b.starts.size() != rows + 1)
      throw std::runtime_error("block " + std::to_string(block) + " has " + std::to_string(b.starts.size() - 1) +
                               " names for " + std::to_string(rows) + " records");
   return b;
}

void NameBlockReader::Load(std::size_t block)
{
   // Read ahead only while the blocks come in order, so random lookups start nothing.
   m_run = (m_block != kNoBlock && block == m_block + 1) ? m_run + 1 : 0;
   const auto rows = [this](std::size_t b) {
      return static_cast<std::size_t>(m_ends[b] - (b == 0 ? 0 : m_ends[b - 1] + 1) + 1);
   };

   auto ahead = m_ahead.find(block);
   if (ahead != m_ahead.end()) {
      m_current = ahead->second.get();
      m_ahead.erase(ahead);
   } else {
      m_current = Decode(Read(block), block, rows(block));
   }
   m_block = block;
   m_firstRow = block == 0 ? 0 : m_ends[block - 1] + 1;

   const std::size_t depth = std::min(m_run, kNameReadAhead);
   for (std::size_t next = block + 1; next < m_ends.size() && next <= block + depth; next++) {
      if (m_ahead.count(next) == 0)
         m_ahead.emplace(next, std::async(std::launch::async, Decode, Read(next), next, rows(next)));
   }
}

SeqBlockReader::SeqBlockReader(RAMFile &file)
   : m_consensusView(file.GetReader().GetView<std::string>(RAMNTupleRecord::kConsensusField)),
     m_posView(file.GetReader().GetView<std::int32_t>(RAMNTupleRecord::kConsensusPosField)),
     m_ends(file.GetMetadata().block_ends)
{
   CheckBlockEnds(m_ends, file.GetReader().GetNEntries());
}

std::string_view SeqBlockReader::Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row)
{
   if (!rec.TestBit(RAMNTupleRecord::kSeqConsensus))
      return rec.seq_lit;
   if (m_block == kNoBlock || row < m_firstRow || row > m_ends[m_block]) {
      const auto it = std::lower_bound(m_ends.begin(), m_ends.end(), static_cast<uint64_t>(row));
      if (it == m_ends.end())
         throw std::runtime_error("record " + std::to_string(row) + " has no consensus");
      m_block = static_cast<std::size_t>(it - m_ends.begin());
      m_firstRow = m_block == 0 ? 0 : m_ends[m_block - 1] + 1;
      m_consensus = m_consensusView(*it);
      m_pos = m_posView(*it);
   }
   Decode(rec, row);
   return m_seq;
}

void SeqBlockReader::Decode(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row)
{
   const auto damaged = [row] {
      return std::runtime_error("the SEQ of record " + std::to_string(row) + " is damaged");
   };
   // One pass over the CIGAR finds where each stretch of SEQ comes from.
   m_pieces.clear();
   std::size_t q = 0;
   std::size_t lit = 0;
   std::int64_t at = static_cast<std::int64_t>(rec.pos) - m_pos;
   const std::vector<uint32_t> &cigar = rec.cigar;
   for (std::size_t i = 0; i < cigar.size();) {
      const int type = CigarType(cigar[i]);
      if (type < 0)
         throw damaged();
      if (type == kReadAndRef) {
         // Long reads have hundreds of =, X and M operations in a row, which
         // read consecutive consensus bases and make one stretch.
         std::size_t n = 0;
         for (; i < cigar.size() && ((kReadAndRefOps >> (cigar[i] & 0xF)) & 1U) != 0; i++)
            n += cigar[i] >> 4;
         if (at < 0)
            throw damaged();
         m_pieces.push_back({q, static_cast<std::size_t>(at), n, false});
         q += n;
         at += static_cast<std::int64_t>(n);
         continue;
      }
      const std::size_t n = cigar[i] >> 4;
      if (type & kRead) {
         if (n > rec.seq_lit.size() - lit)
            throw damaged();
         m_pieces.push_back({q, lit, n, true});
         lit += n;
         q += n;
      }
      if (type & kRef)
         at += static_cast<std::int64_t>(n);
      i++;
   }
   if (lit != rec.seq_lit.size())
      throw damaged();
   m_seq.resize(q);
   for (const Piece &piece : m_pieces) {
      const std::string_view from = piece.literal ? std::string_view(rec.seq_lit) : std::string_view(m_consensus);
      if (piece.from > from.size() || piece.length > from.size() - piece.from)
         throw damaged();
      from.copy(m_seq.data() + piece.to, piece.length, piece.from);
   }

   std::size_t pos = 0;
   for (const uint32_t diff : rec.seq_diff) {
      pos += diff >> kDiffCodeBits;
      if (pos >= m_seq.size())
         throw damaged();
      const char base = kSubstitute[diff & ((1U << kDiffCodeBits) - 1)][static_cast<unsigned char>(m_seq[pos])];
      if (base == '\0')
         throw damaged();
      m_seq[pos] = base;
   }
}
