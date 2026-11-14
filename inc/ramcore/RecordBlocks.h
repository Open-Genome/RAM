#ifndef RAMCORE_RECORDBLOCKS_H
#define RAMCORE_RECORDBLOCKS_H

// Records in blocks. The last record of a block carries the block's names,
// compressed with htscodecs' name tokeniser (tok3, CRAM 3.1's codec), in
// "nameblock", and its lossless qualities in "qualblock", compressed with
// htscodecs' rANS Nx16 as CRAM 3.1 compresses its quality series, or with
// fqzcomp when asked for. It also carries the consensus of the block's aligned
// bases, which the SEQ of its aligned reads is coded against, as CRAM's
// embed_ref=2 slices are, and the bases where the reference its reads' MD tags
// give differs from the consensus. Against that reference MD and NM are
// rebuilt, so a read whose tags come back exactly stores neither. METADATA lists
// the row each block ends on.

#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <ROOT/RNTupleView.hxx>

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <future>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

/// Records per block, as many as a CRAM slice. A region query decodes every
/// block it touches: on NA12878 chr20 about 3 ms for the names, 2 ms for rANS
/// qualities and 60 ms for fqzcomp ones. For names 10,000 gives 60.9 MB,
/// against 62.9 MB with 5,000 and 58.9 MB with 100,000; for fqzcomp it keeps
/// all but a percent of 25,000's gain at less than half the decode time.
constexpr std::size_t kBlockRecords = 10000;
/// Bases, or lossless quality values, per block at most, as CRAM caps a slice
/// at 500 bases per record of kBlockRecords. Without it a block of 10,000
/// PacBio HiFi reads held 150 MB of qualities, all decoded for a region of one
/// base; the writer also holds a whole block's records until it ends.
constexpr std::size_t kBlockBases = 500 * kBlockRecords;
/// Reference bases a block's consensus covers at most, from its leftmost read.
/// A sorted block of short reads covers about 50 kb and one of HiFi reads
/// about 110 kb; the cap bounds the writer's base counts (20 bytes a base) on
/// unsorted or very sparse input, whose reads past it keep SEQ as text.
constexpr std::int64_t kMaxConsensusSpan = 1 << 20;
/// Blocks decoded ahead during a scan. Four keep up with a full dump of fqzcomp
/// blocks, and each fqzcomp decoder holds a 68 MB model, so more only costs
/// memory.
constexpr std::size_t kQualityReadAhead = 4;
/// Blocks of names decoded ahead during a scan. A block decodes in about 3 ms,
/// and a dump writes its records in about 15, so two are plenty.
constexpr std::size_t kNameReadAhead = 2;

/// Holds a block's records until it ends, so that their SEQ can be coded
/// against the block's consensus and the record that ends the block can carry
/// it. Every record's name goes into its block. A quality goes into a block
/// when it is lossless SAM text in '!'..'~'; the rest stay in the record. The
/// block is compressed with fqzcomp when its records' policy has kFqzcomp, and
/// with rANS otherwise.
class BlockWriter {
public:
   using FillFn = std::function<void(ROOT::REntry &)>;

   /// With \p background, a block that ends is coded on another thread while
   /// the next one fills, for a caller that writes on one thread; its rows are
   /// still filled in order on the caller's thread.
   BlockWriter(std::unique_ptr<ROOT::REntry> entry, FillFn fill, std::size_t block_records = kBlockRecords,
               std::size_t block_bases = kBlockBases, bool background = false);
   ~BlockWriter();
   BlockWriter(const BlockWriter &) = delete;
   BlockWriter &operator=(const BlockWriter &) = delete;
   BlockWriter(BlockWriter &&) = delete;
   BlockWriter &operator=(BlockWriter &&) = delete;

   /// The record to fill next. It may hold an earlier record's values.
   RAMNTupleRecord &Record();
   void Add();
   /// Ends the current block and fills every row; needed before the written
   /// rows are flushed.
   void Finish();
   /// Last row of each block finished since the previous call, counted from the
   /// first row after it.
   std::vector<uint64_t> TakeBlockEnds();

private:
   /// A block's records and what is coded from them. Encode() uses nothing
   /// else, so it can run on another thread.
   struct Block {
      std::vector<RAMNTupleRecord> records; ///< Reused: holds at least `rows` records.
      std::size_t rows = 0;
      std::size_t bases = 0;
      std::string names; ///< Each ended by '\0'.
      std::string quals;
      std::vector<uint32_t> lengths;
      std::vector<uint32_t> flags;
      bool fqzcomp = false;                             ///< A record of the block asked for fqzcomp.
      std::vector<std::array<std::uint32_t, 5>> counts; ///< A, C, G, T and others per consensus base.
      std::string consensus;
      std::int32_t consensusPos = 0;
      std::string reference;            ///< The consensus with the bases the reads' MD tags give.
      std::vector<std::uint32_t> patch; ///< Where reference differs from consensus, as kConsensusPatchField.
      std::string md;                   ///< A rebuilt MD, reused between records.
      std::vector<std::uint8_t> qualBlob;
      std::vector<std::uint8_t> nameBlob;
   };

   /// Fills the previous block and encodes this one.
   void EndBlock();
   static void Encode(Block &b);
   /// Sets the records' stored SEQ, the block's consensus and its patches, and
   /// leaves out the MD and NM that can be rebuilt.
   static void CodeSequences(Block &b);
   /// Fills the rows of an encoded block, in order, and empties it.
   void Fill(Block &b);

   std::unique_ptr<ROOT::REntry> m_entry;
   ROOT::RFieldToken m_recordToken;
   std::vector<std::uint8_t> *m_qualBlob = nullptr;
   std::vector<std::uint8_t> *m_nameBlob = nullptr;
   std::string *m_consensus = nullptr;
   std::int32_t *m_consensusPos = nullptr;
   std::vector<std::uint32_t> *m_patch = nullptr;
   FillFn m_fill;
   std::size_t m_blockRecords;
   std::size_t m_blockBases;
   bool m_background;
   Block m_open;   ///< Filled through Record() and Add().
   Block m_closed; ///< Ended; being encoded while m_encoding is valid.
   uint64_t m_rowsSinceTake = 0;
   std::vector<uint64_t> m_blockEnds;
   std::future<void> m_encoding;
};

/// Returns QUAL as SAM text, decoding one block at a time and decoding ahead on
/// other threads during a scan.
class QualityBlockReader {
public:
   explicit QualityBlockReader(RAMFile &file);

   std::string Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row);

private:
   struct Packed {
      std::vector<std::size_t> slots; ///< Per row: index of its quality in the block.
      std::size_t records = 0;
      std::vector<char> bytes;
   };
   struct Block {
      std::vector<std::size_t> slots;
      std::string quals;
      std::vector<int> lengths;
      std::vector<std::size_t> offsets;
   };
   static constexpr std::size_t kNoBlock = static_cast<std::size_t>(-1);

   /// Views are not thread-safe, so only Decode() runs on other threads.
   Packed Read(std::size_t block);
   static Block Decode(Packed packed);
   void Load(std::size_t block);

   ROOT::RNTupleView<uint32_t> m_flagsView;
   std::optional<ROOT::RNTupleView<std::vector<std::uint8_t>>> m_view;
   std::vector<uint64_t> m_ends; ///< Last row of each block.
   std::size_t m_readAhead;
   std::size_t m_block = kNoBlock;
   std::size_t m_run = 0; ///< Consecutive blocks read so far.
   uint64_t m_firstRow = 0;
   Block m_current;
   std::map<std::size_t, std::future<Block>> m_ahead;
};

/// Returns QNAME, decoding one block at a time and decoding ahead on another
/// thread during a scan. Decoding takes about 0.3 us of CPU a name, most of it
/// in the name tokeniser's allocations.
class NameBlockReader {
public:
   explicit NameBlockReader(RAMFile &file);

   /// Valid until the next call.
   std::string_view Get(ROOT::NTupleSize_t row);

private:
   struct Block {
      std::string names;               ///< The block's names, each ended by '\0'.
      std::vector<std::size_t> starts; ///< Where each name starts, and one past the last.
   };
   static constexpr std::size_t kNoBlock = static_cast<std::size_t>(-1);

   /// Views are not thread-safe, so only Decode() runs on other threads.
   std::vector<std::uint8_t> Read(std::size_t block);
   static Block Decode(std::vector<std::uint8_t> packed, std::size_t block, std::size_t rows);
   void Load(std::size_t block);

   ROOT::RNTupleView<std::vector<std::uint8_t>> m_view;
   std::vector<uint64_t> m_ends; ///< Last row of each block.
   std::size_t m_block = kNoBlock;
   std::size_t m_run = 0; ///< Consecutive blocks read so far.
   uint64_t m_firstRow = 0;
   Block m_current;
   std::map<std::size_t, std::future<Block>> m_ahead;
};

/// Returns SEQ, rebuilding the bases its record codes against its block's
/// consensus, and rebuilds the record's MD and NM against the block's reference.
class SeqBlockReader {
public:
   explicit SeqBlockReader(RAMFile &file);

   /// Valid until the next call. Fills the MD and NM of \p rec that the file
   /// left out and clears kMdRebuilt and kNmRebuilt, so call it before
   /// RAMFile::AppendTags().
   std::string_view Get(RAMNTupleRecord &rec, ROOT::NTupleSize_t row);

private:
   static constexpr std::size_t kNoBlock = static_cast<std::size_t>(-1);

   /// A stretch of SEQ copied from the consensus or from seq_lit.
   struct Piece {
      std::size_t to;
      std::size_t from;
      std::size_t length;
      bool literal;
   };

   /// Rebuilds the SEQ of \p rec, row \p row, into m_seq.
   void Decode(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row);

   ROOT::RNTupleView<std::string> m_consensusView;
   ROOT::RNTupleView<std::int32_t> m_posView;
   ROOT::RNTupleView<std::vector<std::uint32_t>> m_patchView;
   std::vector<uint64_t> m_ends; ///< Last row of each block.
   std::size_t m_block = kNoBlock;
   uint64_t m_firstRow = 0;
   std::string m_consensus;
   std::string m_reference;              ///< m_consensus with the block's patches.
   std::vector<std::uint32_t> m_patches; ///< The block's patches, ascending.
   std::int32_t m_pos = 0;
   std::vector<Piece> m_pieces;
   std::string m_seq;
   std::string m_md; ///< A rebuilt MD the record does not keep, when only NM is rebuilt.
};

#endif // RAMCORE_RECORDBLOCKS_H
