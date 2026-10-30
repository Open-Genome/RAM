#include "ramcore/BamtoNTuple.h"
#include "ramcore/RecordBlocks.h"

#include "rntuple/RAMNTupleRecord.h"
#include "rntuple/RAMTags.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <TFile.h>
#include <TStopwatch.h>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdint>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

std::string GetSeq(const bam1_t *b)
{
   const int len = b->core.l_qseq;
   // A BAM record with no sequence is SEQ "*" in SAM, not an empty column.
   if (len <= 0)
      return "*";
   const uint8_t *enc = bam_get_seq(b);
   std::string s(static_cast<std::size_t>(len), 'N');
   for (int i = 0; i < len; ++i)
      s[static_cast<std::size_t>(i)] = seq_nt16_str[bam_seqi(enc, i)];
   return s;
}

std::string GetQual(const bam1_t *b)
{
   const int len = b->core.l_qseq;
   if (len <= 0)
      return "*";
   const uint8_t *q = bam_get_qual(b);
   if (q == nullptr || q[0] == 0xff)
      return "*";
   std::string s(static_cast<std::size_t>(len), '!');
   for (int i = 0; i < len; ++i)
      s[static_cast<std::size_t>(i)] = static_cast<char>(q[i] + 33);
   return s;
}

std::string GetCigar(const bam1_t *b)
{
   const uint32_t n = b->core.n_cigar;
   if (n == 0)
      return "*";
   const uint32_t *cigar = bam_get_cigar(b);
   static constexpr std::array<char, 9> ops = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};
   std::string s;
   s.reserve(n * 4U);
   for (uint32_t i = 0; i < n; ++i) {
      s += std::to_string(bam_cigar_oplen(cigar[i]));
      s += ops.at(static_cast<std::size_t>(bam_cigar_op(cigar[i])));
   }
   return s;
}

// Every aux field as SAM text, formatted by htslib as samtools view prints it.
std::vector<std::string> GetTags(const bam1_t *b)
{
   std::vector<std::string> tags;
   kstring_t text = KS_INITIALIZE;
   const std::unique_ptr<kstring_t, void (*)(kstring_t *)> freeText(&text, [](kstring_t *ks) { ks_free(ks); });
   for (const uint8_t *s = bam_aux_first(b); s != nullptr; s = bam_aux_next(b, s)) {
      ks_clear(&text);
      // The key, type and value around s, as htslib's own bam_aux_get_str() passes them.
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      if (sam_format_aux1(s - 2, *s, s + 1, b->data + b->l_data, &text) == nullptr)
         throw std::runtime_error(std::string("cannot read the tags of ") + bam_get_qname(b));
      tags.emplace_back(text.s, text.l);
   }
   // The walk ends with ENOENT after the last field and EINVAL on a damaged one.
   if (errno == EINVAL)
      throw std::runtime_error(std::string("damaged tags in ") + bam_get_qname(b));
   return tags;
}

void FillRecord(RAMNTupleRecord *rec, const bam1_t *b, const sam_hdr_t *hdr, uint32_t quality_policy,
                RAMNTupleRefs &rnames, RAMNTupleRefs &rnexts, RAMTagWriter &tags)
{
   const char *rname = (b->core.tid < 0) ? "*" : sam_hdr_tid2name(hdr, b->core.tid);

   const char *rnext = "*";
   if (b->core.mtid >= 0)
      rnext = (b->core.mtid == b->core.tid) ? "=" : sam_hdr_tid2name(hdr, b->core.mtid);

   rec->SetBit(quality_policy);
   rec->SetQNAME(bam_get_qname(b));
   rec->SetFLAG(b->core.flag);
   rec->refid = rnames.GetRefId(rname);
   rec->SetPOS(static_cast<int32_t>(b->core.pos + 1));
   rec->SetMAPQ(b->core.qual);
   rec->SetCIGAR(GetCigar(b));
   rec->refnext = rnexts.GetRefId(rnext);
   rec->SetMate(static_cast<int32_t>(b->core.mpos + 1), static_cast<int32_t>(b->core.isize),
                RAMMateOnSameRef(rname, rnext));
   rec->SetSEQ(GetSeq(b));
   rec->SetQUAL(GetQual(b));
   tags.Split(GetTags(b), *rec);
}

} // namespace

void bamtoramntuple(const char *bamfile, const char *treefile, bool /*split*/, bool /*cache*/,
                    int compression_algorithm, uint32_t quality_policy)
{
   TStopwatch stopwatch;
   stopwatch.Start();

   samFile *bamIn = sam_open(bamfile, "r");
   if (bamIn == nullptr) {
      std::cerr << "Cannot open BAM: " << bamfile << "\n";
      return;
   }

   hts_set_threads(bamIn, static_cast<int>(std::thread::hardware_concurrency()));

   sam_hdr_t *hdr = sam_hdr_read(bamIn);
   if (hdr == nullptr) {
      std::cerr << "Cannot read BAM header\n";
      sam_close(bamIn);
      return;
   }

   auto rootFile = std::unique_ptr<TFile>(TFile::Open(treefile, "RECREATE"));
   if (rootFile == nullptr || !rootFile->IsOpen()) {
      std::cerr << "Failed to create: " << treefile << "\n";
      sam_hdr_destroy(hdr);
      sam_close(bamIn);
      return;
   }

   // The state of the file being written: filled in as the records are written,
   // stored in METADATA at the end.
   RAMNTupleRefs rnames;
   RAMNTupleRefs rnexts;
   RAMTagTables tagTables;
   RAMTagWriter tags(tagTables);
   RAMFileMetadata meta;
   RAMCoordinateOrder order;

   const int nRef = sam_hdr_nref(hdr);
   for (int i = 0; i < nRef; ++i)
      rnames.GetRefId(sam_hdr_tid2name(hdr, i));

   auto model = RAMNTupleRecord::MakeModel();
   ROOT::RNTupleWriteOptions opts;
   opts.SetCompression(compression_algorithm);
   opts.SetMaxUnzippedPageSize(64000);

   auto writer = ROOT::RNTupleWriter::Append(std::move(model), "RAM", *rootFile, opts);
   BlockWriter out(writer->GetModel().CreateEntry(), writer->GetModel().CreateEntry(),
                   [&writer](ROOT::REntry &e) { writer->Fill(e); });

   bam1_t *rec = bam_init1();
   std::size_t count = 0;

   while (sam_read1(bamIn, hdr, rec) >= 0) {
      RAMNTupleRecord &record = out.Record();
      FillRecord(&record, rec, hdr, quality_policy, rnames, rnexts, tags);
      meta.max_ref_span = std::max(meta.max_ref_span, record.GetRefSpan());
      order.Note(record.refid, record.pos);
      meta.counts.Note(record.refid, record.flag);
      out.Add();

      ++count;
      if (count % 1000000 == 0)
         std::cerr << "Converted " << count << " records...\n";
   }

   bam_destroy1(rec);
   out.Finish();
   meta.block_ends = out.TakeBlockEnds();
   writer.reset();

   // Region queries can only seek on a sorted file; the file records which it is.
   meta.coordinate_sorted = order.sorted;
   if (!meta.coordinate_sorted)
      std::cerr << bamfile << " is not in coordinate order; region queries will read it in full.\n";

   const char *text = sam_hdr_str(hdr);
   if (text != nullptr) {
      std::istringstream ss(text);
      std::string line;
      while (std::getline(ss, line)) {
         if (line.size() >= 3)
            meta.header_lines.push_back(line);
      }
   }
   meta.rname_refs = rnames.GetRefs();
   meta.rnext_refs = rnexts.GetRefs();
   tagTables.Store(meta);
   RAMNTupleRecord::WriteMetadata(*rootFile, meta);

   rootFile->Close();
   sam_hdr_destroy(hdr);
   sam_close(bamIn);

   std::cout << "\nRAM file created: " << treefile << "\n"
             << "Number of entries: " << count << "\n";
   rnames.Print();
   rnexts.Print();

   stopwatch.Print();
}