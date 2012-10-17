#ifndef AFFINESWOPTIMIZATIONWRAPPER_H_
#define AFFINESWOPTIMIZATIONWRAPPER_H_

#ifdef __cplusplus 
extern "C" {
#endif

#ifdef __cplusplus 
    typedef class AffineSWOptimization tmap_vsw_wrapper_t; 
#else
    typedef struct AffineSWOptimization tmap_vsw_wrapper_t; 
#endif

    tmap_vsw_wrapper_t* 
      tmap_vsw_wrapper_init(int type);

    int32_t
      tmap_vsw_wrapper_process(tmap_vsw_wrapper_t *v,
                               const uint8_t *target, int32_t tlen,
                               const uint8_t *query, int32_t qlen,
                               int32_t mm, int32_t mi, int32_t o, int32_t e, int32_t dir,
                               int32_t qsc, int32_t qec,
                               int32_t *opt, int32_t *target_end, int32_t *query_start, int32_t *n_best);

    void
      tmap_vsw_wrapper_destroy(tmap_vsw_wrapper_t *v);

    int
      tmap_vsw_wrapper_get_max_qlen(tmap_vsw_wrapper_t *v);
    
    int
      tmap_vsw_wrapper_get_max_tlen(tmap_vsw_wrapper_t *v);
#ifdef __cplusplus 
}
#endif

#endif
