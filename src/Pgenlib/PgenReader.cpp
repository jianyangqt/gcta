/* A reader library for plink2 PGEN format
 * This code assembled and revised from plink2.
 * This library removes lots of dependencies of plink2 codes, keep clean and minimal
 * The PGEN format seems to be a draft, thus use the library they provided would be a safe way
 * Use plink2 for QC on genotype in BED or PGEN format.
 * The library don't provide any validity function, only support PGEN generated from plink2
 * Limits:  only support bialleric currently, some bugs in multialleric data
 *
 * Assembled and renewed by Zhili Zheng <zhilizheng@outlook.com>
 * Bug report to Zhili Zheng
 * Copyright (c) 2019 Zhili Zheng
 * Please refer to plink2 for orginal license statement and authorship
 * https://github.com/chrchang/plink-ng
 * Copyright (c) 2005-2019 Shaun Purcell, Christopher Chang.
 *
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "PgenReader.h"

#include <iostream>
#include <limits>

#include "pgenlib_misc.h"
#include "plink2_base.h"

RefcountedWptr *CreateRefcountedWptr(uintptr_t size) {
  RefcountedWptr *rwp = S_CAST(RefcountedWptr *, malloc(sizeof(RefcountedWptr) + size * sizeof(intptr_t)));
  if (!rwp) {
    return nullptr;
  }
  rwp->ref_ct = 1;
  return rwp;
}

void CondReleaseRefcountedWptr(RefcountedWptr **rwpp) {
  RefcountedWptr *rwp = *rwpp;
  if (!rwp) {
    return;
  }
  --rwp->ref_ct;
  if (!rwp->ref_ct) {
    free(rwp);
  }
  *rwpp = nullptr;
}

void stop(const char *prompt) {
  std::cout << prompt << std::endl;
  exit(1);
}

PgenReader::PgenReader() : _info_ptr(nullptr),
                           _allele_idx_offsetsp(nullptr),
                           _nonref_flagsp(nullptr),
                           _state_ptr(nullptr) {
}

void PgenReader::Load(string filename, uint32_t *raw_sample_ct, uint32_t *raw_marker_ct, const vector<uint32_t> &sample_subset_0based) {
  if (_info_ptr) {
    Close();
  }
  _info_ptr = static_cast<plink2::PgenFileInfo *>(malloc(sizeof(plink2::PgenFileInfo)));
  if (!_info_ptr) {
    stop("Out of memory");
  }
  plink2::PreinitPgfi(_info_ptr);
  uint32_t cur_sample_ct = UINT32_MAX;
  if (raw_sample_ct) {
    cur_sample_ct = *raw_sample_ct;
  }
  uint32_t cur_variant_ct = UINT32_MAX;
  if (raw_marker_ct) {
    cur_variant_ct = *raw_marker_ct;
  }
  const char *fname = filename.c_str();
  plink2::PgenHeaderCtrl header_ctrl;
  uintptr_t pgfi_alloc_cacheline_ct;
  char errstr_buf[plink2::kPglErrstrBufBlen];
  if (plink2::PgfiInitPhase1(fname, "", cur_variant_ct, cur_sample_ct, &header_ctrl, _info_ptr, &pgfi_alloc_cacheline_ct, errstr_buf) != plink2::kPglRetSuccess) {
    stop(&(errstr_buf[7]));
  }
  const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
  /*
     if (pvar.isNotNull()) {
     List pvarl = as<List>(pvar);
     if (strcmp_r_c(pvarl[0], "pvar")) {
     stop("pvar is not a pvar object");
     }
     XPtr<class RPvar> rp = as<XPtr<class RPvar>>(pvarl[1]);
     if (rp->GetVariantCt() != raw_variant_ct) {
     stop("pvar and pgen have different variant counts");
     }
     _allele_idx_offsetsp = rp->GetAlleleIdxOffsetsp();
     if (_allele_idx_offsetsp) {
     _info_ptr->allele_idx_offsets = _allele_idx_offsetsp->p;
     }
     _info_ptr->max_allele_ct = rp->GetMaxAlleleCt();
     } else {
     */
  if (header_ctrl & 0x30) {
    // no need to zero-initialize this
    _allele_idx_offsetsp = CreateRefcountedWptr(raw_variant_ct + 1);
    _info_ptr->allele_idx_offsets = _allele_idx_offsetsp->p;
    // _info_ptr->max_allele_ct updated by PgfiInitPhase2() in this case
  }
  _info_ptr->max_allele_ct = 2;
  //}

  if ((header_ctrl & 0xc0) == 0xc0) {
    // todo: load this in pvar, to enable consistency check.  we use a
    // (manually implemented) shared_ptr in preparation for this.
    const uintptr_t raw_variant_ctl = plink2::DivUp(raw_variant_ct, plink2::kBitsPerWord);
    // no need to zero-initialize this
    _nonref_flagsp = CreateRefcountedWptr(raw_variant_ctl + 1);
    _info_ptr->nonref_flags = _nonref_flagsp->p;
  }
  const uint32_t file_sample_ct = _info_ptr->raw_sample_ct;
  _raw_sample_ct = file_sample_ct;
  unsigned char *pgfi_alloc = nullptr;
  if (plink2::cachealigned_malloc(pgfi_alloc_cacheline_ct * plink2::kCacheline, &pgfi_alloc)) {
    stop("Out of memory");
  }
  uint32_t max_vrec_width;
  uintptr_t pgr_alloc_cacheline_ct;
  if (PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, raw_variant_ct, &max_vrec_width, _info_ptr, pgfi_alloc, &pgr_alloc_cacheline_ct, errstr_buf)) {
    if (pgfi_alloc && (!_info_ptr->vrtypes)) {
      plink2::aligned_free(pgfi_alloc);
    }
    stop(&(errstr_buf[7]));
  }
  if ((!_allele_idx_offsetsp) && (_info_ptr->gflags & 4)) {
    // Note that it's safe to be ignorant of multiallelic variants when
    // phase and dosage info aren't present; GetAlleleCt() then always returns
    // 2 when that isn't actually true, and all ALTs are treated as if they
    // were ALT1, but otherwise everything works properly.
    stop("Multiallelic variants and phase/dosage info simultaneously present; pvar required in this case");
  }
  _state_ptr = static_cast<plink2::PgenReader *>(malloc(sizeof(plink2::PgenReader)));
  if (!_state_ptr) {
    stop("Out of memory");
  }
  plink2::PreinitPgr(_state_ptr);
  plink2::PgrSetFreadBuf(nullptr, _state_ptr);
  const uintptr_t pgr_alloc_main_byte_ct = pgr_alloc_cacheline_ct * plink2::kCacheline;
  const uintptr_t sample_subset_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerVec) * plink2::kBytesPerVec;
  const uintptr_t cumulative_popcounts_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerWord * plink2::kInt32PerVec) * plink2::kBytesPerVec;
  const uintptr_t genovec_byte_ct = plink2::DivUp(file_sample_ct, plink2::kNypsPerVec) * plink2::kBytesPerVec;
  _genovecbuf_uptr_size = genovec_byte_ct;
  _genovec_rawbuf_uptr_size = genovec_byte_ct;
  const uintptr_t ac_byte_ct = plink2::RoundUpPow2(file_sample_ct * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
  const uintptr_t ac2_byte_ct = plink2::RoundUpPow2(file_sample_ct * 2 * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
  uintptr_t multiallelic_hc_byte_ct = 0;
  if (_info_ptr->max_allele_ct != 2) {
    multiallelic_hc_byte_ct = 2 * sample_subset_byte_ct + ac_byte_ct + ac2_byte_ct;
  }
  const uintptr_t dosage_main_byte_ct = plink2::DivUp(file_sample_ct, (2 * plink2::kInt32PerVec)) * plink2::kBytesPerVec;
  unsigned char *pgr_alloc;
  if (plink2::cachealigned_malloc(pgr_alloc_main_byte_ct + (2 * plink2::kPglNypTransposeBatch + 5) * sample_subset_byte_ct + cumulative_popcounts_byte_ct + (1 + plink2::kPglNypTransposeBatch) * genovec_byte_ct + multiallelic_hc_byte_ct + dosage_main_byte_ct + plink2::kPglBitTransposeBufbytes + 4 * (plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8), &pgr_alloc)) {
    stop("Out of memory");
  }
  plink2::PglErr reterr = PgrInit(fname, max_vrec_width, _info_ptr, _state_ptr, pgr_alloc);
  if (reterr != plink2::kPglRetSuccess) {
    if (!plink2::PgrGetFreadBuf(_state_ptr)) {
      plink2::aligned_free(pgr_alloc);
    }
    sprintf(errstr_buf, "PgrInit() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
  unsigned char *pgr_alloc_iter = &(pgr_alloc[pgr_alloc_main_byte_ct]);
  _subset_include_vec = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _subset_include_interleaved_vec = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);

#ifdef USE_AVX2
  _subset_include_interleaved_vec[-3] = 0;
  _subset_include_interleaved_vec[-2] = 0;
#endif
  _subset_include_interleaved_vec[-1] = 0;

  _subset_cumulative_popcounts = reinterpret_cast<uint32_t *>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[cumulative_popcounts_byte_ct]);
  _pgv.genovec = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[genovec_byte_ct]);
  if (multiallelic_hc_byte_ct) {
    _pgv.patch_01_set = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.patch_01_vals = reinterpret_cast<plink2::AlleleCode *>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[ac_byte_ct]);
    _pgv.patch_10_set = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.patch_10_vals = reinterpret_cast<plink2::AlleleCode *>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[ac2_byte_ct]);
  } else {
    _pgv.patch_01_set = nullptr;
    _pgv.patch_01_vals = nullptr;
    _pgv.patch_10_set = nullptr;
    _pgv.patch_10_vals = nullptr;
  }
  _pgv.phasepresent = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _pgv.phaseinfo = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _pgv.dosage_present = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _pgv.dosage_main = reinterpret_cast<uint16_t *>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[dosage_main_byte_ct]);
  // std::cout << "phase info size: " << dosage_main_byte_ct << ", dosage_present: " << sample_subset_byte_ct << std::endl;;
  if (sample_subset_0based.size() > 0) {
    SetSampleSubsetInternal(sample_subset_0based);
  } else {
    _subset_size = file_sample_ct;
  }
  _transpose_batch_buf = reinterpret_cast<plink2::VecW *>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglBitTransposeBufbytes]);
  _multivar_vmaj_geno_buf = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * genovec_byte_ct]);
  _multivar_vmaj_phasepresent_buf = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * sample_subset_byte_ct]);
  _multivar_vmaj_phaseinfo_buf = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * sample_subset_byte_ct]);
  _multivar_smaj_geno_batch_buf = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 4]);
  _multivar_smaj_phaseinfo_batch_buf = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8]);
  _multivar_smaj_phasepresent_batch_buf = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
  // pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8]);
  // for dosage
  genovec_size64 = (GetGenoBufPtrSize(_subset_size) + 63) / 64 * 64;
  dosage_present_size64 = (GetDosagePresentSize(_subset_size) + 63) / 64 * 64;
  dosage_main_size64 = (GetDosageMainSize(_subset_size) + 63) / 64 * 64;
}

uint32_t PgenReader::GetRawSampleCt() const {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  return _info_ptr->raw_sample_ct;
}

uint32_t PgenReader::GetSubsetSize() const {
  return _subset_size;
}

uint32_t PgenReader::GetVariantCt() const {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  return _info_ptr->raw_variant_ct;
}

uint32_t PgenReader::GetAlleleCt(uint32_t variant_idx) const {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if (variant_idx >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  if (!_allele_idx_offsetsp) {
    return 2;
  }
  const uintptr_t *allele_idx_offsets = _allele_idx_offsetsp->p;
  return allele_idx_offsets[variant_idx + 1] - allele_idx_offsets[variant_idx];
}

uint32_t PgenReader::GetMaxAlleleCt() const {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  return _info_ptr->max_allele_ct;
}

uintptr_t PgenReader::GetGenoBufSizeUptr() const {
  return _genovecbuf_uptr_size;
}

uintptr_t PgenReader::GetGenoBufFullSizeUptr() const {
  return _genovec_rawbuf_uptr_size;
}

bool PgenReader::HardcallPhasePresent() const {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  return ((_info_ptr->gflags & plink2::kfPgenGlobalHardcallPhasePresent) != 0);
}
// kfPgenGlobalDosagePresent

// true: have phase, false: don't have phase
bool PgenReader::IsPhasePresent() const {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  return ((_info_ptr->gflags & (plink2::kfPgenGlobalDosagePhasePresent | plink2::kfPgenGlobalHardcallPhasePresent)) != 0);
}

bool PgenReader::HasMultiAllelic() const {
  return ((_info_ptr->gflags & plink2::kfPgenGlobalMultiallelicHardcallFound) != 0);
}

static const int32_t kGenoRInt32Quads[1024] ALIGNV16 = QUAD_TABLE256(0, 1, 2, -9);

void PgenReader::ReadIntHardcalls(vector<int32_t> &buf, int variant_idx, int allele_idx) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if (static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  if (buf.size() != _subset_size) {
    char errstr_buf[256];
    sprintf(errstr_buf, "buf has wrong length (%" PRIdPTR "; %u expected)", buf.size(), _subset_size);
    stop(errstr_buf);
  }
  plink2::PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(_subset_cumulative_popcounts, _state_ptr, &pssi);
  plink2::PglErr reterr = plink2::PgrGet1(_subset_include_vec, pssi, _subset_size, variant_idx, allele_idx, _state_ptr, _pgv.genovec);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "PgrGet1() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
  plink2::GenoarrLookup256x4bx4(_pgv.genovec, kGenoRInt32Quads, _subset_size, buf.data());
}

static const double kGenoRDoublePairs[32] ALIGNV16 = PAIR_TABLE16(0.0, 1.0, 2.0, std::numeric_limits<double>::quiet_NaN());

void PgenReader::ReadRawHard(uintptr_t *buf, int variant_idx, int allele_idx) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if (static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  /*
  if (buf.size() != _genovecbuf_uptr_size) {
      stop("buffer size not enough");
  }
  */
  plink2::PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(_subset_cumulative_popcounts, _state_ptr, &pssi);
  plink2::PglErr reterr = plink2::PgrGet1(_subset_include_vec, pssi, _subset_size, variant_idx, allele_idx, _state_ptr, buf);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "PgrGet1() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
}

void PgenReader::ReadRawFullHard(uintptr_t *buf, int variant_idx, int allele_idx) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if (static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  /*
  if (buf.size() != _genovec_rawbuf_uptr_size) {
      stop("buffer size not correct");
  }
  */
  // cout << "start read" << std::endl;
  plink2::PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(nullptr, _state_ptr, &pssi);
  plink2::PglErr reterr = plink2::PgrGet1(NULL, pssi, _info_ptr->raw_sample_ct, variant_idx, allele_idx, _state_ptr, buf);
  /* // not much faster
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  plink2::PgenReaderMain * pgrp = &(_state_ptr->GET_PRIVATE_m());
  plink2::PglErr reterr = plink2::ReadRawGenovec(false, variant_idx, pgrp, &fread_ptr, &fread_end, buf.data());
  */

  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "PgrGet1() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
}

void PgenReader::ReadRawFullHard(uintptr_t *buf, int variant_idx) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if (static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  plink2::PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(nullptr, _state_ptr, &pssi);
  plink2::PglErr reterr = plink2::PgrGet(NULL, pssi, _info_ptr->raw_sample_ct, variant_idx,  _state_ptr, buf);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "PgrGet() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
}

bool PgenReader::CountHardFreqMissExt(uintptr_t *buf, const uintptr_t *subset_iter_vec, uint32_t rawSampleSize, uint32_t keepSize, SNPInfo *snpinfo, bool f_std) {
  std::array<uint32_t, 4> genocounts;
  plink2::ZeroTrailingNyps(rawSampleSize, buf);

  if (rawSampleSize == keepSize) {
    uint32_t sample_ct = keepSize;
    // const uint32_t sample_ct_remainder = sample_ct % plink2::kBitsPerWordD2;
    plink2::GenoarrCountFreqsUnsafe(buf, sample_ct, genocounts);
    /*
    plink2::GenoarrCountFreqsUnsafe(buf, sample_ct - sample_ct_remainder, genocounts);
    if (sample_ct_remainder) {
        uintptr_t cur_geno_word = plink2::bzhi(buf[sample_ct / plink2::kBitsPerWordD2], 2 * sample_ct_remainder);
        const uintptr_t cur_geno_word_high = plink2::kMask5555 & (cur_geno_word >> 1);
        const uint32_t even_ct = plink2::Popcount01Word(cur_geno_word & plink2::kMask5555);
        const uint32_t odd_ct = plink2::Popcount01Word(cur_geno_word_high);
        const uint32_t bothset_ct = plink2::Popcount01Word(cur_geno_word & cur_geno_word_high);
        genocounts[0] += sample_ct_remainder + bothset_ct - even_ct - odd_ct;
        genocounts[1] += even_ct - bothset_ct;
        genocounts[2] += odd_ct - bothset_ct;
        genocounts[3] += bothset_ct;
    }
    */
  } else {
    const uint32_t raw_sample_ct = rawSampleSize;
    plink2::GenoarrCountSubsetFreqs(buf, subset_iter_vec, raw_sample_ct, keepSize, genocounts);
  }
  snpinfo->N = keepSize - genocounts[3];
  snpinfo->AlCount = 2 * snpinfo->N;
  uint32_t dGenoCounts2 = 2 * genocounts[2];
  uint32_t dosage = genocounts[1] + dGenoCounts2;
  snpinfo->af = 1.0 * dosage / snpinfo->AlCount;
  snpinfo->mean = 2.0 * snpinfo->af;
  if (f_std) snpinfo->std = (dosage + dGenoCounts2 - snpinfo->mean * dosage) / (snpinfo->N - 1);
  snpinfo->nMissRate = 1.0 * snpinfo->N / keepSize;
  return true;
  /*
     std::cout << "00: " << genocounts[0] << "\n";
     std::cout << "01: " << genocounts[1] << "\n";
     std::cout << "10: " << genocounts[2] << "\n";
     std::cout << "11: " << genocounts[3] << "\n";
     std::cout << "tt: " << genocounts[0] + genocounts[1]
     + genocounts[2] + genocounts[3] << "\n";
     */
}
#include <csignal>
bool PgenReader::CountHardDosage(uintptr_t *buf, uint16_t *dosage_buf, const uintptr_t *dosage_present, const vector<uint32_t> *maskp, uint32_t sampleCT, uint32_t dosageCT, SNPInfo *snpinfo, string &err) {
  if (sampleCT != dosageCT) {
    std::raise(SIGINT);
    std::cout << sampleCT << ", dosageCT: " << dosageCT << std::endl;
    err = "GCTA can't support mix hard-call and dosage currently.";
    return false;
  }
  uint64_t sum = 0, sumsq = 0;
  for (int i = 0; i < dosageCT; i++) {
    uint16_t curDos = dosage_buf[i];
    sum += curDos;
    sumsq += curDos * curDos;
  }
  uint32_t maskSize = 0;
  if (maskp != NULL) {
    const vector<uint32_t> masks = *maskp;
    maskSize = masks.size();
    for (int i = 0; i < maskSize; i++) {
      uint16_t curDos = dosage_buf[masks[i]];
      sum -= curDos / 2;
      sumsq -= 3 * curDos / 4;
    }
  }
  snpinfo->nMissRate = 1;
  snpinfo->N = sampleCT;
  snpinfo->AlCount = 2 * sampleCT - maskSize;
  snpinfo->af = 0.00006103515625 * sum / snpinfo->AlCount;
  snpinfo->std = 0.00006103515625 * 0.00006103515625 * (1.0 * sumsq - 1.0 * sum * sum / sampleCT) / (sampleCT - 1);

  return true;
}

bool PgenReader::CountHardFreqMissExtX(uintptr_t *buf, const uintptr_t *subset_iter_vec, const uintptr_t *subset_iter_vec2,
                                       uint32_t rawSampleSize, uint32_t keepSize, uint32_t keepSize2, SNPInfo *snpinfo, string &errmsg, bool dosageComp, bool f_std) {
  std::array<uint32_t, 4> genocounts;
  std::array<uint32_t, 4> genocounts2;
  plink2::ZeroTrailingNyps(rawSampleSize, buf);
  if (rawSampleSize == keepSize) {
    uint32_t sample_ct = keepSize;
    plink2::GenoarrCountFreqsUnsafe(buf, sample_ct, genocounts);
    if (keepSize == keepSize2) {
      // plink2::GenoarrCountFreqsUnsafe(buf, sample_ct, genocounts);
      genocounts2 = genocounts;
    } else {
      plink2::GenoarrCountSubsetFreqs(buf, subset_iter_vec2, rawSampleSize, keepSize2, genocounts2);
    }
  } else {
    plink2::GenoarrCountSubsetFreqs(buf, subset_iter_vec, rawSampleSize, keepSize, genocounts);
    if (keepSize2 != 0) {
      plink2::GenoarrCountSubsetFreqs(buf, subset_iter_vec2, rawSampleSize, keepSize2, genocounts2);
    } else {
      genocounts2 = {0, 0, 0, 0};
    }
  }
  bool retVal = true;
  if (genocounts2[1] != 0) {
    errmsg = "het";
    retVal = false;
  }
  snpinfo->N = keepSize - genocounts[3];
  snpinfo->AlCount = 2 * snpinfo->N - genocounts2[0] - genocounts2[1] - genocounts2[2];
  snpinfo->nMissRate = 1.0 * snpinfo->N / keepSize;
  // double dosage = genocounts[1] + 2 * genocounts[2] - genocounts2[2] - 0.5 * genocounts2[1];
  int dGenoCount1Homo = 2 * genocounts[2];
  double dosage = genocounts[1] + dGenoCount1Homo;
  double dosage_half = dosage - genocounts2[2] - 0.5 * genocounts2[1];

  snpinfo->af = dosage_half / snpinfo->AlCount;
  if (dosageComp) {
    snpinfo->mean = dosage / snpinfo->N;
    if (f_std) snpinfo->std = (dosage + dGenoCount1Homo - dosage * snpinfo->mean) / (snpinfo->N - 1);
  } else {
    snpinfo->mean = dosage_half / snpinfo->N;
    if (f_std) snpinfo->std = (dosage + dGenoCount1Homo - 0.75 * genocounts2[1] - 3 * genocounts2[2] - dosage_half * snpinfo->mean) / (snpinfo->N - 1);
  }
  return retVal;
}

void PgenReader::CountHardFreqMiss(uintptr_t *buf, SNPInfo *snpinfo) {
  PgenReader::CountHardFreqMissExt(buf, _subset_include_interleaved_vec, _info_ptr->raw_sample_ct, _subset_size, snpinfo, true);
}

void PgenReader::ExtractGeno(const uintptr_t *in, uintptr_t *out) {
  PgenReader::ExtractGenoExt(in, _subset_include_vec, _raw_sample_ct, _subset_size, out);
}

void PgenReader::ExtractGenoExt(const uintptr_t *in, const uintptr_t *subsets, uint32_t rawSampleSize, uint32_t keepSize, uintptr_t *out) {
  plink2::CopyNyparrNonemptySubset(in, subsets, rawSampleSize, keepSize, out);
}

void PgenReader::ExtractDoubleExt(uintptr_t *in, const uintptr_t *subsets, uint32_t rawSampleSize, uint32_t keepSize, const double *gtable, double *gOut, uintptr_t *missOut) {
  uintptr_t *bufptr = NULL;
  bool newBuf = false;
  if (rawSampleSize == keepSize) {
    bufptr = in;
  } else {
    newBuf = true;
    bufptr = new uintptr_t[GetGenoBufPtrSize(keepSize)];
    ExtractGenoExt(in, subsets, rawSampleSize, keepSize, bufptr);
  }
  plink2::GenoarrLookup16x8bx2(bufptr, gtable, keepSize, gOut);
  if (missOut != NULL) {
    plink2::GenoarrToMissingnessUnsafe(bufptr, keepSize, missOut);
  }
  if (newBuf) {
    delete[] bufptr;
  }
}

void PgenReader::ReadHardcalls(vector<double> &buf, int variant_idx, int allele_idx) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if (static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  if (buf.size() != _subset_size) {
    char errstr_buf[256];
    sprintf(errstr_buf, "buf has wrong length (%" PRIdPTR "; %u expected)", buf.size(), _subset_size);
    stop(errstr_buf);
  }
  plink2::PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(_subset_cumulative_popcounts, _state_ptr, &pssi);
  plink2::PglErr reterr = plink2::PgrGet1(_subset_include_vec, pssi, _subset_size, variant_idx, allele_idx, _state_ptr, _pgv.genovec);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "PgrGet1() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
  plink2::GenoarrLookup16x8bx2(_pgv.genovec, kGenoRDoublePairs, _subset_size, buf.data());
}

void PgenReader::Read(vector<double> &buf, int variant_idx, int allele_idx) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if (static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  if (buf.size() != _subset_size) {
    char errstr_buf[256];
    sprintf(errstr_buf, "buf has wrong length (%" PRIdPTR "; %u expected)", buf.size(), _subset_size);
    stop(errstr_buf);
  }
  uint32_t dosage_ct;
  plink2::PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(_subset_cumulative_popcounts, _state_ptr, &pssi);
  plink2::PglErr reterr = plink2::PgrGet1D(_subset_include_vec, pssi, _subset_size, variant_idx, allele_idx, _state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "PgrGet1D() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
  plink2::Dosage16ToDoubles(kGenoRDoublePairs, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, buf.data());
}

void PgenReader::ReadDosage(uintptr_t *buf, int variant_idx, int allele_idx) {
  uintptr_t *genovec = buf;
  uintptr_t *dosage_present = (buf + genovec_size64);
  // uintptr_t *dosage_present = buf;
  // uintptr_t *genovec = buf + dosage_present_size64;
  uint16_t *dosage_main = reinterpret_cast<uint16_t *>(buf + genovec_size64 + dosage_present_size64);
  uint32_t *dosage_ct_ptr = reinterpret_cast<uint32_t *>(buf + genovec_size64 + dosage_present_size64 + dosage_main_size64);
  // std::cout << "genovec_size64: " << genovec_size64 << ", dosage present: " << dosage_present_size64 << ", dosage main: " << dosage_main_size64 << std::endl;
  // plink2::PglErr reterr = PgrGet1D(_subset_include_vec, _subset_cumulative_popcounts, _subset_size, variant_idx, allele_idx, _state_ptr, genovec, dosage_present, dosage_main, dosage_ct_ptr);
  uint32_t dosage_ct;
  plink2::PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(_subset_cumulative_popcounts, _state_ptr, &pssi);
  plink2::PglErr reterr = plink2::PgrGet1D(_subset_include_vec, pssi, _subset_size, variant_idx, allele_idx, _state_ptr, genovec, dosage_present, dosage_main, dosage_ct_ptr);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "PgrGet1D() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
  // plink2::Dosage16ToDoubles(kGenoRDoublePairs, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, buf.data());
}

static const uint64_t kGenoToRIntcodeDPairs[32] ALIGNV16 = PAIR_TABLE16(0, 0x100000000LLU, 0x100000001LLU, 0x8000000080000000LLU);
static const int32_t kGenoToLogicalPhaseQuads[1024] ALIGNV16 = QUAD_TABLE256(1, 0, 1, -9);
/*
   void PgenReader::ReadAlleles(IntegerMatrix acbuf, Nullable<LogicalVector> phasepresent_buf, int variant_idx) {
   if (!_info_ptr) {
   stop("pgen is closed");
   }
   if ((acbuf.nrow() != 2) || (acbuf.ncol() != static_cast<int>(_subset_size))) {
   char errstr_buf[256];
   sprintf(errstr_buf, "acbuf has wrong size (%dx%d; 2x%u expected)", acbuf.nrow(), acbuf.ncol(), _subset_size);
   stop(errstr_buf);
   }
   ReadAllelesPhasedInternal(variant_idx);
   plink2::GenoarrToAlleleCodes(kGenoToRIntcodeDPairs, _pgv.genovec, _subset_size, &acbuf[0]);
   const uintptr_t* allele_idx_offsets = _info_ptr->allele_idx_offsets;
   uint32_t cur_allele_ct = 2;
   if (allele_idx_offsets) {
   cur_allele_ct = allele_idx_offsets[variant_idx + 1] - allele_idx_offsets[variant_idx];
   if (cur_allele_ct != 2) {
   stop("multiallelic support under development");
   }
   }
   const uintptr_t* phasepresent = _pgv.phasepresent;
   const uintptr_t* phaseinfo = _pgv.phaseinfo;
   const uint32_t phasepresent_ct = _pgv.phasepresent_ct;
   uintptr_t sample_uidx_base = 0;
   uintptr_t cur_bits = phasepresent[0];
   if (!phasepresent_buf.isNotNull()) {
   if (cur_allele_ct == 2) {
   uint64_t* allele_codes_alias64 = R_CAST(uint64_t*, &acbuf[0]);
   for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
   const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
   if (plink2::IsSet(phaseinfo, sample_uidx)) {
// 1|0
allele_codes_alias64[sample_uidx] = 1;
}
}
} else {
int32_t* allele_codes = &acbuf[0];
for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
if (plink2::IsSet(phaseinfo, sample_uidx)) {
const int32_t tmpval = allele_codes[2 * sample_uidx];
allele_codes[2 * sample_uidx] = allele_codes[2 * sample_uidx + 1];
allele_codes[2 * sample_uidx + 1] = tmpval;
}
}
}
return;
}
// Unfortunately, we can't use GenoarrPhasedToAlleleCodes directly, since
// it's written for Python 1-byte bools instead of R 4-byte logical values.
// (probable todo: allow the no-phasepresent_buf part to be called
// separately)
//
// 0, 2 -> automatically phased.  3 -> NA_LOGICAL.
// 1 -> assume unphased; then change to phased as necessary when iterating
//      over phasepresent.
int32_t* phasepresent_wbuf = &(as<LogicalVector>(phasepresent_buf)[0]);
plink2::GenoarrLookup256x4bx4(_pgv.genovec, kGenoToLogicalPhaseQuads, _subset_size, phasepresent_wbuf);
if (cur_allele_ct == 2) {
uint64_t* allele_codes_alias64 = R_CAST(uint64_t*, &acbuf[0]);
for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
phasepresent_wbuf[sample_uidx] = 1;
if (plink2::IsSet(phaseinfo, sample_uidx)) {
allele_codes_alias64[sample_uidx] = 1;
}
}
} else {
int32_t* allele_codes = &acbuf[0];
for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
phasepresent_wbuf[sample_uidx] = 1;
if (plink2::IsSet(phaseinfo, sample_uidx)) {
    const int32_t tmpval = allele_codes[2 * sample_uidx];
    allele_codes[2 * sample_uidx] = allele_codes[2 * sample_uidx + 1];
    allele_codes[2 * sample_uidx + 1] = tmpval;
}
}
}
}

static const double kGenoToRNumcodePairs[8] ALIGNV16 = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, NA_REAL, NA_REAL};

void PgenReader::ReadAllelesNumeric(NumericMatrix acbuf, Nullable<LogicalVector> phasepresent_buf, int variant_idx) {
    if (!_info_ptr) {
        stop("pgen is closed");
    }
    if ((acbuf.nrow() != 2) || (acbuf.ncol() != static_cast<int>(_subset_size))) {
        char errstr_buf[256];
        sprintf(errstr_buf, "acbuf has wrong size (%dx%d; 2x%u expected)", acbuf.nrow(), acbuf.ncol(), _subset_size);
        stop(errstr_buf);
    }
    ReadAllelesPhasedInternal(variant_idx);
    double* allele_codes = &acbuf[0];
    plink2::GenoarrLookup4x16b(_pgv.genovec, kGenoToRNumcodePairs, _subset_size, allele_codes);
    const uintptr_t* allele_idx_offsets = _info_ptr->allele_idx_offsets;
    uint32_t cur_allele_ct = 2;
    if (allele_idx_offsets) {
        cur_allele_ct = allele_idx_offsets[variant_idx + 1] - allele_idx_offsets[variant_idx];
        if (cur_allele_ct != 2) {
            stop("multiallelic support under development");
        }
    }
    const uintptr_t* phasepresent = _pgv.phasepresent;
    const uintptr_t* phaseinfo = _pgv.phaseinfo;
    const uint32_t phasepresent_ct = _pgv.phasepresent_ct;
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = phasepresent[0];
    if (!phasepresent_buf.isNotNull()) {
        if (cur_allele_ct == 2) {
            for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
                const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
                if (plink2::IsSet(phaseinfo, sample_uidx)) {
                    // 1|0
                    allele_codes[2 * sample_uidx] = 1.0;
                    allele_codes[2 * sample_uidx + 1] = 0.0;
                }
            }
        } else {
            for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
                const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
                if (plink2::IsSet(phaseinfo, sample_uidx)) {
                    const double tmpval = allele_codes[2 * sample_uidx];
                    allele_codes[2 * sample_uidx] = allele_codes[2 * sample_uidx + 1];
                    allele_codes[2 * sample_uidx + 1] = tmpval;
                }
            }
        }
        return;
    }
    int32_t* phasepresent_wbuf = &(as<LogicalVector>(phasepresent_buf)[0]);
    plink2::GenoarrLookup256x4bx4(_pgv.genovec, kGenoToLogicalPhaseQuads, _subset_size, phasepresent_wbuf);
    if (cur_allele_ct == 2) {
        for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
            const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
            phasepresent_wbuf[sample_uidx] = 1;
            if (plink2::IsSet(phaseinfo, sample_uidx)) {
                allele_codes[2 * sample_uidx] = 1.0;
                allele_codes[2 * sample_uidx + 1] = 0.0;
            }
        }
    } else {
        for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
            const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
            phasepresent_wbuf[sample_uidx] = 1;
            if (plink2::IsSet(phaseinfo, sample_uidx)) {
                const double tmpval = allele_codes[2 * sample_uidx];
                allele_codes[2 * sample_uidx] = allele_codes[2 * sample_uidx + 1];
                allele_codes[2 * sample_uidx + 1] = tmpval;
            }
        }
    }
}

void PgenReader::ReadIntList(IntegerMatrix buf, IntegerVector variant_subset) {
    if (!_info_ptr) {
        stop("pgen is closed");
    }
    // assume that buf has the correct dimensions
    const uintptr_t vsubset_size = variant_subset.size();
    const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
    int32_t* buf_iter = &buf[0];
    for (uintptr_t col_idx = 0; col_idx != vsubset_size; ++col_idx) {
        uint32_t variant_idx = variant_subset[col_idx] - 1;
        if (static_cast<uint32_t>(variant_idx) >= raw_variant_ct) {
            char errstr_buf[256];
            sprintf(errstr_buf, "variant_subset element out of range (%d; must be 1..%u)", variant_idx + 1, raw_variant_ct);
            stop(errstr_buf);
        }
        plink2::PglErr reterr = PgrGet(_subset_include_vec, _subset_cumulative_popcounts, _subset_size, variant_idx, _state_ptr, _pgv.genovec);
        if (reterr != plink2::kPglRetSuccess) {
            char errstr_buf[256];
            sprintf(errstr_buf, "PgrGet() error %d", static_cast<int>(reterr));
            stop(errstr_buf);
        }
        plink2::GenoarrLookup256x4bx4(_pgv.genovec, kGenoRInt32Quads, _subset_size, buf_iter);
        buf_iter = &(buf_iter[_subset_size]);
    }
}

void PgenReader::ReadList(NumericMatrix buf, IntegerVector variant_subset, bool meanimpute) {
    if (!_info_ptr) {
        stop("pgen is closed");
    }
    // assume that buf has the correct dimensions
    const uintptr_t vsubset_size = variant_subset.size();
    const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
    double* buf_iter = &buf[0];
    for (uintptr_t col_idx = 0; col_idx != vsubset_size; ++col_idx) {
        uint32_t variant_idx = variant_subset[col_idx] - 1;
        if (static_cast<uint32_t>(variant_idx) >= raw_variant_ct) {
            char errstr_buf[256];
            sprintf(errstr_buf, "variant_subset element out of range (%d; must be 1..%u)", variant_idx + 1, raw_variant_ct);
            stop(errstr_buf);
        }
        uint32_t dosage_ct;
        plink2::PglErr reterr = PgrGetD(_subset_include_vec, _subset_cumulative_popcounts, _subset_size, variant_idx, _state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
        if (reterr != plink2::kPglRetSuccess) {
            char errstr_buf[256];
            sprintf(errstr_buf, "PgrGetD() error %d", static_cast<int>(reterr));
            stop(errstr_buf);
        }
        if (!meanimpute) {
            plink2::Dosage16ToDoubles(kGenoRDoublePairs, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, buf_iter);
        } else {
            plink2::ZeroTrailingNyps(_subset_size, _pgv.genovec);
            if (plink2::Dosage16ToDoublesMeanimpute(_pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, buf_iter)) {
                char errstr_buf[256];
                sprintf(errstr_buf, "variant %d has only missing dosages", variant_idx + 1);
                stop(errstr_buf);
            }
        }
        buf_iter = &(buf_iter[_subset_size]);
    }
}

void PgenReader::FillVariantScores(NumericVector result, NumericVector weights, Nullable<IntegerVector> variant_subset) {
    if (!_info_ptr) {
        stop("pgen is closed");
    }
    if (weights.size() != _subset_size) {
        char errstr_buf[256];
        sprintf(errstr_buf, "weights.size()=%td doesn't match pgen sample-subset size=%d", weights.size(), _subset_size);
        stop(errstr_buf);
    }
    const int raw_variant_ct = _info_ptr->raw_variant_ct;
    const int* variant_idx_ints = nullptr;
    uintptr_t variant_ct = raw_variant_ct;
    if (variant_subset.isNotNull()) {
        IntegerVector vs = as<IntegerVector>(variant_subset);
        variant_idx_ints = &(vs[0]);
        variant_ct = vs.size();
    }
    for (uintptr_t ulii = 0; ulii != variant_ct; ++ulii) {
        int variant_idx = ulii;
        if (variant_idx_ints) {
            variant_idx = variant_idx_ints[ulii] - 1;
            if ((variant_idx < 0) || (variant_idx >= raw_variant_ct)) {
                char errstr_buf[256];
                sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, raw_variant_ct);
                stop(errstr_buf);
            }
        }
        uint32_t dosage_ct;
        plink2::PglErr reterr = plink2::PgrGetD(_subset_include_vec, _subset_cumulative_popcounts, _subset_size, variant_idx, _state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
        if (reterr != plink2::kPglRetSuccess) {
            char errstr_buf[256];
            sprintf(errstr_buf, "PgrGetD() error %d", static_cast<int>(reterr));
            stop(errstr_buf);
        }
        plink2::ZeroTrailingNyps(_subset_size, _pgv.genovec);
        const double* wts = &(weights[0]);
        result[ulii] = plink2::LinearCombinationMeanimpute(wts, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct);
    }
}
*/

void PgenReader::Close() {
  // don't bother propagating file close errors for now
  if (_info_ptr) {
    CondReleaseRefcountedWptr(&_allele_idx_offsetsp);
    CondReleaseRefcountedWptr(&_nonref_flagsp);
    if (_info_ptr->vrtypes) {
      plink2::aligned_free(_info_ptr->vrtypes);
    }
    plink2::PglErr reterr = plink2::kPglRetSuccess;
    plink2::CleanupPgfi(_info_ptr, &reterr);
    free(_info_ptr);
    _info_ptr = nullptr;
  }
  if (_state_ptr) {
    plink2::PglErr reterr = plink2::kPglRetSuccess;
    plink2::CleanupPgr(_state_ptr, &reterr);
    if (PgrGetFreadBuf(_state_ptr)) {
      plink2::aligned_free(PgrGetFreadBuf(_state_ptr));
    }
    free(_state_ptr);
    _state_ptr = nullptr;
  }
  _subset_size = 0;
}

void PgenReader::SetSampleSubsets(const vector<uint32_t> &sample_subset_0based, uint32_t rawSampleSize, uintptr_t *subset_incl_vec, uintptr_t *subset_inter_vec) {
  const uint32_t raw_sample_ct = rawSampleSize;
  // plink2::kBitsPerVec = 128
  // plink2::kWordsPerVec = 2
  const uint32_t raw_sample_ctv = plink2::DivUp(raw_sample_ct, plink2::kBitsPerVec);
  const uint32_t raw_sample_ctaw = raw_sample_ctv * plink2::kWordsPerVec;
  uintptr_t *sample_include = subset_incl_vec;
  plink2::ZeroWArr(raw_sample_ctaw, sample_include);
  const uint32_t subset_size = sample_subset_0based.size();
  if (subset_size == 0) {
    return;
    // stop("Empty sample_subset is not currently permitted");
  }
  uint32_t sample_uidx = sample_subset_0based[0];
  uint32_t idx = 0;
  uint32_t next_uidx;
  while (1) {
    if (sample_uidx >= raw_sample_ct) {
      char errstr_buf[256];
      sprintf(errstr_buf, "sample number out of range (%d; must be 1..%u)", static_cast<int>(sample_uidx + 1), raw_sample_ct);
      stop(errstr_buf);
    }
    plink2::SetBit(sample_uidx, sample_include);
    if (++idx == subset_size) {
      break;
    }
    next_uidx = sample_subset_0based[idx];

    // prohibit this since it implies that the caller expects genotypes to be
    // returned in a different order
    if (next_uidx <= sample_uidx) {
      std::cout << next_uidx << "\t" << sample_uidx << "\t" << idx << "\n";
      stop("sample_subset is not in strictly increasing order");
    }
    sample_uidx = next_uidx;
  }
  plink2::FillInterleavedMaskVec(sample_include, raw_sample_ctv, subset_inter_vec);
}

void PgenReader::SetSampleSubsetInternal(const vector<uint32_t> &sample_subset_0based) {
  const uint32_t raw_sample_ct = _info_ptr->raw_sample_ct;
  PgenReader::SetSampleSubsets(sample_subset_0based, raw_sample_ct, _subset_include_vec, _subset_include_interleaved_vec);

  const uint32_t raw_sample_ctl = plink2::DivUp(raw_sample_ct, plink2::kBitsPerWord);
  plink2::FillCumulativePopcounts(_subset_include_vec, raw_sample_ctl, _subset_cumulative_popcounts);
  _subset_size = sample_subset_0based.size();
  _genovecbuf_uptr_size = plink2::DivUp(_subset_size, plink2::kNypsPerVec) * plink2::kBytesPerVec;
}

void PgenReader::ReadAllelesPhasedInternal(int variant_idx) {
  if (static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  plink2::PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(_subset_cumulative_popcounts, _state_ptr, &pssi);
  plink2::PglErr reterr = plink2::PgrGetMP(_subset_include_vec, pssi, _subset_size, variant_idx, _state_ptr, &_pgv);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "PgrGetMP() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
}
int PgenReader::GetGenoBufPtrSize(uint32_t sample_ct) {
  // plink2::kNypsPerVec = 64
  // plink2::kBytesPerVec = 16
  return plink2::DivUp(sample_ct, plink2::kNypsPerVec) * plink2::kBytesPerVec;
}

int PgenReader::GetSubsetMaskSize(uint32_t sample_ct) {
  const uint32_t raw_sample_ct = sample_ct;
  const uint32_t raw_sample_ctv = plink2::DivUp(raw_sample_ct, plink2::kBitsPerVec);  // plink2::kBitsPerVec = 128
  const uint32_t raw_sample_ctaw = raw_sample_ctv * plink2::kWordsPerVec;             // plink2::kWordsPerVec = 2
  return raw_sample_ctaw;
}

int PgenReader::GetDosagePresentSize(uint32_t sample_ct) {
  const int perSize = sizeof(uintptr_t);  // 8
  // plink2::kBitsPerVec = 128
  // plink2::kBytesPerVec = 16
  return plink2::DivUp(sample_ct, plink2::kBitsPerVec) * plink2::kBytesPerVec / perSize;
}

int PgenReader::GetDosageMainSize(uint32_t sample_ct) {
  const int perSize = sizeof(uintptr_t);  // 8
  // plink2::kInt32PerVec = 4
  // plink2::kBytesPerVec = 16
  return plink2::DivUp(sample_ct, (2 * plink2::kInt32PerVec)) * plink2::kBytesPerVec / perSize;
}

PgenReader::~PgenReader() {
  Close();
}
// plink2::Pack11ToHalfword
