struct FABDriverBase
{
    const double *in_ptr[c_NUM];
    const int *in_lo;
    const int *in_hi;
    int stride[3];

    double *m_out_ptr[c_NUM];
    const int *m_out_lo;
    const int *m_out_hi;
    int m_out_stride[3];

    template <class data_t>
    void local_vars(int idx, vars_t<data_t>& out);

    template <class data_t>
    void diff1(int idx, int stride, vars_t<data_t>& out);

    template <class data_t>
    void diff2(int idx, int stride, vars_t<data_t>& out);

    template <class data_t>
    void mixed_diff2(int idx, int stride1, int stride2, vars_t<data_t>& out);

    template <class data_t>
    void advection(int idx, const tensor<1, data_t>& vec, vars_t<data_t>& out);

    template <class data_t>
    void dissipation(int idx, vars_t<data_t>& out);
};

template <class compute_t>
class FABDriver : public FABDriverBase
{
public:
    FABDriverContext m_ctx;
    compute_t& m_compute;

    template <typename... param_types>
    FABDriver(param_types... params);

    void execute(const FArrayBox& in, FArrayBox& out);
};
