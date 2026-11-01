#include "ramcore/RecordBlocks.h"
#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <TFile.h>

#include <gtest/gtest.h>

#include <cstdint>
#include <cstdio>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Read {
   uint16_t flag;
   std::string qual;
   std::string name = "r";
};

const char *const kFile = "record_blocks_test.root";
// The policies samtoramntuple passes for lossless qualities, by default and
// with -fqzcomp.
const uint32_t kRansPolicy = RAMNTupleRecord::kPhred33;
const uint32_t kFqzcompPolicy = RAMNTupleRecord::kPhred33 | RAMNTupleRecord::kFqzcomp;
// The first byte of a quality block, by codec.
const uint8_t kRansCodec = 1;
const uint8_t kFqzcompCodec = 2;

using Damage = std::function<void(std::vector<std::uint8_t> &)>;

// Writes the reads; \p ends, when given, replaces the block ends, and \p damage
// changes every quality block, as in a damaged file.
void Write(const std::vector<Read> &reads, std::size_t block, uint32_t policy, const std::vector<uint64_t> &ends = {},
           const Damage &damage = nullptr, std::size_t bases = kBlockBases)
{
   std::unique_ptr<TFile> file(TFile::Open(kFile, "RECREATE"));
   RAMFileMetadata meta;
   meta.rname_refs = {"chr1"};
   {
      auto writer = ROOT::RNTupleWriter::Append(RAMNTupleRecord::MakeModel(), "RAM", *file);
      BlockWriter out(
         writer->GetModel().CreateEntry(),
         [&writer, &damage](ROOT::REntry &e) {
            auto &blob = *e.GetPtr<std::vector<std::uint8_t>>(RAMNTupleRecord::kQualBlockField);
            if (damage && !blob.empty())
               damage(blob);
            writer->Fill(e);
         },
         block, bases);
      for (const auto &read : reads) {
         RAMNTupleRecord &rec = out.Record();
         rec.SetBit(policy);
         rec.SetQNAME(read.name);
         rec.SetFLAG(read.flag);
         rec.refid = 0;
         rec.SetQUAL(read.qual);
         out.Add();
      }
      out.Finish();
      meta.block_ends = ends.empty() ? out.TakeBlockEnds() : ends;
   }
   RAMNTupleRecord::WriteMetadata(*file, meta);
   file->Close();
}

void ExpectQualities(const std::vector<Read> &reads, const std::vector<std::size_t> &order)
{
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   auto view = file->GetReader().GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*file);
   for (const std::size_t row : order)
      EXPECT_EQ(quals.Get(view(row), row), reads[row].qual) << "row " << row;
}

// The first byte of the quality block that row \p end carries, or 0 for none.
uint8_t CodecAt(uint64_t end)
{
   auto file = RAMFile::Open(kFile);
   EXPECT_NE(file, nullptr);
   if (!file)
      return 0;
   auto blocks = file->GetReader().GetView<std::vector<std::uint8_t>>(RAMNTupleRecord::kQualBlockField);
   const std::vector<std::uint8_t> &blob = blocks(end);
   return blob.empty() ? 0 : blob.front();
}

class RecordBlocksTest : public ::testing::Test {
protected:
   // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
   void TearDown() override { std::remove(kFile); }
};

TEST_F(RecordBlocksTest, EveryQualityComesBackInAnyOrder)
{
   const std::vector<Read> reads = {
      {0, "IIIIHHHGGF"},
      {0x10, "#####ABCDE"},
      {0x80, "!~!~!~"},
      {0, "*"}, // reuses row 1's record object, so its in-block bit must be cleared
      {0x90, "ABCDEFGHIJKLMNO"},
      {0, "AB C"}, // not SAM quality
      {0, "FFFFFFFFFF"},
      {0x10, "5"},
      {0, "*"},
      {0, "JJJJJJJJJJJJ"},
   };
   for (const uint32_t policy : {kRansPolicy, kFqzcompPolicy}) {
      Write(reads, /*block=*/3, policy);
      const auto written = RAMFile::Open(kFile);
      ASSERT_NE(written, nullptr);
      EXPECT_EQ(written->GetMetadata().block_ends, (std::vector<uint64_t>{2, 5, 8, 9}));

      std::vector<std::size_t> forward;
      std::vector<std::size_t> backward;
      for (std::size_t i = 0; i < reads.size(); i++) {
         forward.push_back(i);
         backward.push_back(reads.size() - 1 - i);
      }
      ExpectQualities(reads, forward);
      ExpectQualities(reads, backward);

      auto file = RAMFile::Open(kFile);
      auto view = file->GetReader().GetView<RAMNTupleRecord>("record");
      for (std::size_t row = 0; row < reads.size(); row++) {
         const bool inBlock = reads[row].qual != "*" && reads[row].qual != "AB C";
         EXPECT_EQ(view(row).TestBit(RAMNTupleRecord::kQualInBlock), inBlock) << "row " << row;
         EXPECT_FALSE(view(row).TestBit(RAMNTupleRecord::kFqzcomp)) << "a request, not kept in the record";
         if (inBlock)
            EXPECT_TRUE(view(row).qual.empty()) << "row " << row;
      }
   }
}

// A block also ends once it holds that many quality values, as a CRAM slice
// of long reads ends on bases, and the qualities come back across the ends.
TEST_F(RecordBlocksTest, BlocksOfLongReadsEndOnBases)
{
   const std::vector<Read> reads = {{0, std::string(20, 'I')}, {0x10, std::string(20, '#')}, {0, "*"},
                                    {0, std::string(5, 'F')},  {0, std::string(30, '5')},    {0, "AB"}};
   Write(reads, /*block=*/100, kRansPolicy, /*ends=*/{}, /*damage=*/nullptr, /*bases=*/25);
   const auto written = RAMFile::Open(kFile);
   ASSERT_NE(written, nullptr);
   EXPECT_EQ(written->GetMetadata().block_ends, (std::vector<uint64_t>{1, 4, 5}));
   ExpectQualities(reads, {5, 0, 3, 1, 4, 2});
}

// Lossless qualities go into blocks by default, compressed with rANS; with
// -fqzcomp the same blocks are fqzcomp's.
TEST_F(RecordBlocksTest, EachBlockSaysItsCodec)
{
   const std::vector<Read> reads = {{0, "IIII"}, {0x10, "#+5?"}, {0, "FFFF"}};
   Write(reads, /*block=*/2, kRansPolicy);
   EXPECT_EQ(CodecAt(1), kRansCodec);
   EXPECT_EQ(CodecAt(2), kRansCodec);
   Write(reads, /*block=*/2, kFqzcompPolicy);
   EXPECT_EQ(CodecAt(1), kFqzcompCodec);
   EXPECT_EQ(CodecAt(2), kFqzcompCodec);
}

TEST_F(RecordBlocksTest, ABlockWithoutQualitiesNeedsNoData)
{
   const std::vector<Read> reads = {{0, "*"}, {0, "*"}, {0, "*"}, {0, "HHHH"}};
   for (const uint32_t policy : {kRansPolicy, kFqzcompPolicy}) {
      Write(reads, /*block=*/3, policy);
      EXPECT_EQ(CodecAt(2), 0) << "no qualities, no block";
      ExpectQualities(reads, {3, 0, 1, 2});
   }
}

TEST_F(RecordBlocksTest, LossyPoliciesStayInTheRecord)
{
   for (const uint32_t policy :
        std::vector<uint32_t>{RAMNTupleRecord::kDrop, RAMNTupleRecord::kDrop | RAMNTupleRecord::kFqzcomp}) {
      Write({{0, "IIII"}}, /*block=*/3, policy);
      auto file = RAMFile::Open(kFile);
      ASSERT_NE(file, nullptr);
      auto view = file->GetReader().GetView<RAMNTupleRecord>("record");
      QualityBlockReader quals(*file);
      EXPECT_FALSE(view(0).TestBit(RAMNTupleRecord::kQualInBlock));
      EXPECT_EQ(quals.Get(view(0), 0), "*");
   }
}

// Binned qualities from Q30 up encode to printable bytes ('F', Q37, becomes
// byte 37, '%'), which must not be taken for SAM quality and put in a block.
TEST_F(RecordBlocksTest, IlluminaBinnedQualitiesStayInTheRecord)
{
   Write({{0, "FFFFFFFF"}, {0, "#+5?FFII"}}, /*block=*/3, RAMNTupleRecord::kIlluminaBinning);
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   auto view = file->GetReader().GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*file);
   for (std::size_t row = 0; row < 2; row++) {
      EXPECT_FALSE(view(row).TestBit(RAMNTupleRecord::kQualInBlock)) << "row " << row;
      EXPECT_EQ(quals.Get(view(row), row), view(row).GetQUAL()) << "row " << row;
   }
   EXPECT_EQ(quals.Get(view(0), 0), "FFFFFFFF");
}

// A damaged file is reported, not read past: block ends that fall back, and
// an end beyond the last row.
TEST_F(RecordBlocksTest, DamagedBlockEndsAreReported)
{
   const std::vector<Read> reads = {{0, "IIII"}, {0, "HHHH"}, {0, "GGGG"}};
   for (const auto &ends : {std::vector<uint64_t>{1, 0}, std::vector<uint64_t>{1, 3}}) {
      Write(reads, /*block=*/2, kRansPolicy, ends);
      auto file = RAMFile::Open(kFile);
      ASSERT_NE(file, nullptr);
      EXPECT_THROW(QualityBlockReader quals(*file), std::runtime_error) << ends[0] << "," << ends[1];
      EXPECT_THROW(NameBlockReader names(*file), std::runtime_error) << ends[0] << "," << ends[1];
   }
}

// A damaged quality block is reported when it is read: an unknown codec, a
// block cut short, a lengths part longer than the block, and streams that do
// not decode.
TEST_F(RecordBlocksTest, DamagedQualityBlocksAreReported)
{
   const std::vector<Read> reads = {{0, "IIIIHHHGGF"}, {0x10, "#####ABCDE"}};
   const std::vector<std::pair<uint32_t, Damage>> damages = {
      {kRansPolicy, [](std::vector<std::uint8_t> &b) { b.at(0) = 9; }},
      {kRansPolicy, [](std::vector<std::uint8_t> &b) { b.resize(1); }},
      {kRansPolicy, [](std::vector<std::uint8_t> &b) { b.resize(4); }},
      {kRansPolicy, [](std::vector<std::uint8_t> &b) { b.at(4) = 0x7f; }},
      {kRansPolicy, [](std::vector<std::uint8_t> &b) { b.resize(b.size() - 1); }},
      {kRansPolicy, [](std::vector<std::uint8_t> &b) { b.at(1) = 1; }},
      {kFqzcompPolicy, [](std::vector<std::uint8_t> &b) { b.resize(3); }},
   };
   for (std::size_t i = 0; i < damages.size(); i++) {
      Write(reads, /*block=*/2, damages[i].first, /*ends=*/{}, damages[i].second);
      auto file = RAMFile::Open(kFile);
      ASSERT_NE(file, nullptr);
      auto view = file->GetReader().GetView<RAMNTupleRecord>("record");
      QualityBlockReader quals(*file);
      EXPECT_THROW(quals.Get(view(0), 0), std::runtime_error) << "damage " << i;
   }
}

// Names go into blocks whatever the quality policy, and come back in any order:
// a repeated name, "*", an empty one, the longest SAM allows.
TEST_F(RecordBlocksTest, EveryNameComesBackInAnyOrder)
{
   const std::vector<Read> reads = {
      {0, "IIII", "A00217:76:HFLT3DSXX:4:1428:30599:8797"},
      {0, "IIII", "A00217:76:HFLT3DSXX:4:1428:30599:8797"},
      {0, "IIII", "*"},
      {0, "IIII", ""},
      {0, "IIII", std::string(254, 'q')},
      {0, "IIII", "SRR123.4"},
      {0, "IIII", "r!#$%&'()*+,-./0;<=>?@[]^_`{|}~"},
   };
   Write(reads, /*block=*/3, RAMNTupleRecord::kPhred33);
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   NameBlockReader names(*file);
   for (const std::size_t row : std::vector<std::size_t>{0, 1, 2, 3, 4, 5, 6, 6, 3, 0, 5, 1})
      EXPECT_EQ(names.Get(row), reads[row].name) << "row " << row;
}

TEST_F(RecordBlocksTest, NamesAreNotAColumn)
{
   Write({{0, "IIII", "a"}}, /*block=*/3, RAMNTupleRecord::kPhred33);
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   const auto &desc = file->GetReader().GetDescriptor();
   EXPECT_EQ(desc.FindFieldId("qname", desc.FindFieldId("record")), ROOT::kInvalidDescriptorId);
}

// Block ends that put more rows in a block than it has names for.
TEST_F(RecordBlocksTest, ABlockWithTooFewNamesIsReported)
{
   const std::vector<Read> reads = {{0, "IIII", "a"}, {0, "IIII", "b"}, {0, "IIII", "c"}};
   Write(reads, /*block=*/2, RAMNTupleRecord::kPhred33, /*ends=*/{2});
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   NameBlockReader names(*file);
   EXPECT_THROW(names.Get(0), std::runtime_error);
}

struct Aligned {
   int32_t refid;
   int32_t pos; ///< 0-based
   uint16_t flag;
   std::string cigar;
   std::string seq;
   bool coded = false; ///< Expected to be coded against the consensus, in one block.
};

// Writes the reads through BlockWriter, \p block records to a block, each
// named by its row.
void WriteSeqs(const std::vector<Aligned> &reads, std::size_t block, bool background = false)
{
   std::unique_ptr<TFile> file(TFile::Open(kFile, "RECREATE"));
   RAMFileMetadata meta;
   meta.rname_refs = {"chr1", "chr2"};
   {
      auto writer = ROOT::RNTupleWriter::Append(RAMNTupleRecord::MakeModel(), "RAM", *file);
      BlockWriter out(
         writer->GetModel().CreateEntry(), [&writer](ROOT::REntry &e) { writer->Fill(e); }, block, kBlockBases,
         background);
      for (std::size_t i = 0; i < reads.size(); i++) {
         const Aligned &read = reads[i];
         RAMNTupleRecord &rec = out.Record();
         rec.SetCompressionMode(RAMNTupleRecord::kPhred33);
         rec.SetQNAME(std::to_string(i));
         rec.SetFLAG(read.flag);
         rec.refid = read.refid;
         rec.pos = read.pos;
         rec.SetCIGAR(read.cigar);
         rec.SetSEQ(read.seq);
         rec.SetQUAL("*");
         out.Add();
      }
      out.Finish();
      meta.block_ends = out.TakeBlockEnds();
   }
   RAMNTupleRecord::WriteMetadata(*file, meta);
   file->Close();
}

const std::vector<Aligned> kAligned = {
   {0, 100, 0, "10M", "ACGTACGTAC", true},
   {0, 100, 0, "10M", "ACGTACGTAC", true},
   {0, 100, 0, "10M", "ACGTACGTAC", true},
   // Substitutions of A, C, G and T, and an N and an IUPAC base.
   {0, 100, 0, "10M", "TCATACGNAY", true},
   // Clips, an insertion, a deletion and padding.
   {0, 102, 0, "2H3S4M2I3M1D2M1P1M5H", "GGGGTACTTGTACAC", true},
   {0, 100, 0x10, "4=1X5=", "ACGTTCGTAC", true},
   {0, 100, 0, "2M50N3M", "ACTTT", true},
   // Runs past the others; at 110 it ties G with the clipped read's C, and
   // ties go to A, C, G and T in that order.
   {0, 108, 0, "6M", "ACGGCA", true},
   {0, 100, 0, "10M", "*"},
   {0, 100, 0x4, "*", "ACGT"},
   {-1, -1, 0x4, "*", "TTTT"},
   // The CIGAR has five bases, SEQ seven.
   {0, 100, 0, "5M", "ACGTACG"},
   {1, 100, 0, "4M", "GGGG"},
   {0, 100 + static_cast<int32_t>(kMaxConsensusSpan), 0, "4M", "CCCC"},
   {0, 101, 0, "4M", ""},
};

void ExpectSeqs(const std::vector<Aligned> &reads, bool checkCoded)
{
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   auto view = file->GetReader().GetView<RAMNTupleRecord>("record");
   SeqBlockReader seqs(*file);
   std::vector<std::size_t> order;
   for (std::size_t i = 0; i < reads.size(); i++)
      order.push_back(i);
   for (std::size_t i = reads.size(); i-- > 0;)
      order.push_back(i);
   order.push_back(3);
   for (const std::size_t row : order) {
      const RAMNTupleRecord &rec = view(row);
      EXPECT_EQ(seqs.Get(rec, row), reads[row].seq) << "row " << row;
      if (checkCoded) {
         EXPECT_EQ(rec.TestBit(RAMNTupleRecord::kSeqConsensus), reads[row].coded) << "row " << row;
      }
   }
}

TEST_F(RecordBlocksTest, SeqComesBackFromTheConsensus)
{
   WriteSeqs(kAligned, kBlockRecords);
   ExpectSeqs(kAligned, /*checkCoded=*/true);
}

TEST_F(RecordBlocksTest, SeqComesBackFromSmallBlocks)
{
   for (const std::size_t block : {1U, 2U, 3U, 4U, 5U}) {
      SCOPED_TRACE("block of " + std::to_string(block));
      WriteSeqs(kAligned, block);
      ExpectSeqs(kAligned, /*checkCoded=*/false);
   }
}

// Coded on another thread, the blocks still fill their rows in order.
TEST_F(RecordBlocksTest, BlocksCodedInTheBackgroundKeepTheirOrder)
{
   for (const std::size_t block : {1U, 2U, 3U, 4U, 5U, 100U}) {
      SCOPED_TRACE("block of " + std::to_string(block));
      WriteSeqs(kAligned, block, /*background=*/true);
      ExpectSeqs(kAligned, /*checkCoded=*/block == 100);
      auto file = RAMFile::Open(kFile);
      ASSERT_NE(file, nullptr);
      NameBlockReader names(*file);
      for (std::size_t row = 0; row < kAligned.size(); row++)
         EXPECT_EQ(names.Get(row), std::to_string(row));
   }
}

TEST_F(RecordBlocksTest, TheConsensusTakesTheCommonestBase)
{
   WriteSeqs(kAligned, kBlockRecords);
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   auto consensus = file->GetReader().GetView<std::string>(RAMNTupleRecord::kConsensusField);
   auto pos = file->GetReader().GetView<std::int32_t>(RAMNTupleRecord::kConsensusPosField);
   EXPECT_EQ(pos(kAligned.size() - 1), 100);
   // 114 to 151 lie inside the 50N skip.
   EXPECT_EQ(consensus(kAligned.size() - 1), "ACGTACGTACCACA" + std::string(38, 'N') + "TTT");
   EXPECT_EQ(consensus(0), "");
}

// One record per block, with \p consensus at 0 and a 4M CIGAR, and each of
// \p records as its seq_diff.
void WriteCoded(const std::vector<std::vector<uint32_t>> &records, const std::string &consensus)
{
   std::unique_ptr<TFile> file(TFile::Open(kFile, "RECREATE"));
   RAMFileMetadata meta;
   meta.rname_refs = {"chr1"};
   {
      auto writer = ROOT::RNTupleWriter::Append(RAMNTupleRecord::MakeModel(), "RAM", *file);
      auto entry = writer->GetModel().CreateEntry();
      auto rec = entry->GetPtr<RAMNTupleRecord>("record");
      *entry->GetPtr<std::string>(RAMNTupleRecord::kConsensusField) = consensus;
      for (const auto &diffs : records) {
         rec->SetCompressionMode(RAMNTupleRecord::kSeqConsensus);
         rec->refid = 0;
         rec->pos = 0;
         rec->SetCIGAR("4M");
         rec->seq_diff = diffs;
         writer->Fill(*entry);
         meta.block_ends.push_back(meta.block_ends.size());
      }
   }
   RAMNTupleRecord::WriteMetadata(*file, meta);
   file->Close();
}

TEST_F(RecordBlocksTest, DamagedSeqIsReported)
{
   // A step and a code: (step << 5) | code.
   WriteCoded({{1U << 5}, {4U << 5}, {31U}, {(1U << 5) | 3U, 0U}, {0U}}, "ACGT");
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   auto view = file->GetReader().GetView<RAMNTupleRecord>("record");
   SeqBlockReader seqs(*file);
   EXPECT_EQ(seqs.Get(view(0), 0), "AGGT");
   for (const ROOT::NTupleSize_t row : {1U, 2U})
      EXPECT_THROW(seqs.Get(view(row), row), std::runtime_error) << "row " << row;
   // Code 3 is '=' (kBases[0]), and the base it replaces is no longer A, C, G
   // or T for the next code 0 at the same place.
   EXPECT_THROW(seqs.Get(view(3), 3), std::runtime_error);
   EXPECT_EQ(seqs.Get(view(4), 4), "CCGT");

   WriteCoded({{}}, "ACG");
   auto shortFile = RAMFile::Open(kFile);
   ASSERT_NE(shortFile, nullptr);
   auto shortView = shortFile->GetReader().GetView<RAMNTupleRecord>("record");
   SeqBlockReader shortSeqs(*shortFile);
   EXPECT_THROW(shortSeqs.Get(shortView(0), 0), std::runtime_error);
}

} // namespace
