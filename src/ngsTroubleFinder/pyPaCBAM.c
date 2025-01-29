#include "htslib/hts.h"
#include "htslib/khash.h"
#include "htslib/sam.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

const uint32_t qualityFlags =
    BAM_FDUP | BAM_FQCFAIL | BAM_FUNMAP | BAM_FSECONDARY;
KHASH_SET_INIT_STR(str)

typedef struct {
  htsFile *htf;
  bam_hdr_t *header;
  hts_idx_t *index;

} ngsFile;

int getIndex(bam1_t *read, uint32_t position, uint32_t quality) {
  uint32_t referenceOffset = read->core.pos;
  uint32_t queryOffset = 0;
  uint32_t *cigar = bam_get_cigar(read);

  // The reads starts after the position so it's not mapped to the reference
  if (referenceOffset > position)
    return -1;

  for (uint32_t i = 0; i < read->core.n_cigar; i++) {
    uint32_t element = cigar[i];
    uint32_t operation = element & BAM_CIGAR_MASK;
    uint32_t operationLength = element >> BAM_CIGAR_SHIFT;
    uint32_t consumesQue = (bam_cigar_type(operation) & 1) != 0;
    uint32_t consumesRef = (bam_cigar_type(operation) & 2) != 0;

    if (consumesQue)
      queryOffset += operationLength;

    if (consumesRef)
      referenceOffset += operationLength;

    // We reached the variant position
    if (referenceOffset > position) {
      // Check that we haven't reached it through an indel
      if (consumesQue & consumesRef) {
        uint32_t index = queryOffset - referenceOffset + position;

        // Check that the quality is at least as high as required
        if (bam_get_qual(read)[index] >= quality)
          return index;
        else
          return -1;
      }

      else
        return -1;
    }
  }
  return -1;
}

int inferHaplotype(ngsFile *ngsFile, char *chrom, int positionV1,
                   int positionV2, char ref1, char alt1, char ref2, char alt2,
                   int baseQuality, int readQuality, int *haplotypes) {
  // haplotypes is an int[5]
  // haplotypes represent RefRef, RefAlt, AltRef, AltAlt, Other

  int tid = bam_name2id(ngsFile->header, chrom);
  if (tid == -1) {
    fprintf(stderr,
            "[PyPaCBAM]: Chromosome not found.\n[PyPaCBAM]: Check if vcf "
            "chromosome and bam alignment chromosome name match.\n");
    return 1;
  }

  bam1_t *read = bam_init1();
  hts_itr_t *itr =
      sam_itr_queryi(ngsFile->index, tid, positionV1, positionV2 + 1);
  if (itr == NULL) {
    fprintf(stderr, "[PyPaCBAM]: Cannot retrieve genomic region.\n");
    return 1;
  }

  khash_t(str) * readsAlreadySeen;
  readsAlreadySeen = kh_init(str);

  while (sam_itr_next(ngsFile->htf, itr, read) >= 0) {
    if ((read->core.qual >= readQuality) &
        ((read->core.flag & qualityFlags) == 0) & (read->core.tid >= 0)) {
      char *readName = bam_get_qname(read);
      khint_t k = kh_get(str, readsAlreadySeen, readName);
      int isMissing = (k == kh_end(readsAlreadySeen));
      if (isMissing) {
        int indexV1 = getIndex(read, positionV1, baseQuality);
        int indexV2 = getIndex(read, positionV2, baseQuality);
        if ((indexV1 > -1) & (indexV2 > -1)) {
          // counts the bases based on the following string from sam spec file
          char baseV1 =
              "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(read), indexV1)];
          char baseV2 =
              "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(read), indexV2)];

          if ((baseV1 == ref1) & (baseV2 == ref2))
            haplotypes[0]++;
          else if ((baseV1 == ref1) & (baseV2 == alt2))
            haplotypes[1]++;
          else if ((baseV1 == alt1) & (baseV2 == ref2))
            haplotypes[2]++;
          else if ((baseV1 == alt1) & (baseV2 == alt2))
            haplotypes[3]++;
          else
            haplotypes[4]++;

          int absent;
          k = kh_put(str, readsAlreadySeen, readName, &absent);
          if (absent)
            kh_key(readsAlreadySeen, k) = strdup(readName);
        }
      }
    }
  }

  for (khint_t k = 0; k < kh_end(readsAlreadySeen); k++) {
    if (kh_exist(readsAlreadySeen, k))
      free((char *)kh_key(readsAlreadySeen, k));
  }

  kh_destroy(str, readsAlreadySeen);
  bam_destroy1(read);
  hts_itr_destroy(itr);

  return 0;
}

int pileup(ngsFile *ngsFile, char *chrom, int position, int baseQuality,
           int readQuality, int *pileupResults) {
  int counter[16] = {0};
  int tid = bam_name2id(ngsFile->header, chrom);
  if (tid == -1) {
    fprintf(stderr,
            "[PyPaCBAM]: Chromosome not found.\n[PyPaCBAM]: Check if vcf "
            "chromosome and bam alignment chromosome name match.\n");
    return 1;
  }

  bam1_t *read = bam_init1();
  hts_itr_t *itr = sam_itr_queryi(ngsFile->index, tid, position, position + 1);
  if (itr == NULL) {
    fprintf(stderr, "[PyPaCBAM]: Cannot retrieve genomic region.");
    return 1;
  }

  khash_t(str) * readsAlreadySeen;
  readsAlreadySeen = kh_init(str);

  while (sam_itr_next(ngsFile->htf, itr, read) >= 0) {
    if ((read->core.qual >= readQuality) &
        ((read->core.flag & qualityFlags) == 0) & (read->core.tid >= 0)) {
      char *readName = bam_get_qname(read);
      khint_t k = kh_get(str, readsAlreadySeen, readName);
      int isMissing = (k == kh_end(readsAlreadySeen));
      if (isMissing) {
        int index = getIndex(read, position, baseQuality);
        if (index > -1) {
          // counts the bases based on the following string from sam spec file
          //"=ACMGRSVTWYHKDBN"[base];
          counter[bam_seqi(bam_get_seq(read), index)]++;
          int absent;
          k = kh_put(str, readsAlreadySeen, readName, &absent);
          if (absent)
            kh_key(readsAlreadySeen, k) = strdup(readName);
        }
      }
    }
  }
  for (khint_t k = 0; k < kh_end(readsAlreadySeen); k++) {
    if (kh_exist(readsAlreadySeen, k))
      free((char *)kh_key(readsAlreadySeen, k));
  }

  kh_destroy(str, readsAlreadySeen);
  bam_destroy1(read);
  hts_itr_destroy(itr);

  // Repack the array by getting only the base values
  pileupResults[0] = counter[1];
  pileupResults[1] = counter[2];
  pileupResults[2] = counter[4];
  pileupResults[3] = counter[8];
  return 0;
}

int openHtsFile(char *fileName, char *fileIndex, char *referenceGenome,
                ngsFile *f) {
  f->htf = hts_open(fileName, "rb");
  if (f->htf == NULL)
    return 1;

  if (referenceGenome != NULL)
    hts_set_opt(f->htf, CRAM_OPT_REFERENCE, referenceGenome);

  f->header = sam_hdr_read(f->htf);
  if (f->header == NULL)
    return 1;

  f->index = sam_index_load(f->htf, fileIndex);
  if (f->index == NULL)
    return 1;

  return 0;
}

void closeHtsFile(ngsFile *f) {
  bam_hdr_destroy(f->header);
  hts_idx_destroy(f->index);
  hts_close(f->htf);
}
