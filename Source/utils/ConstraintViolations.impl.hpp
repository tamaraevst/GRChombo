/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CONSTRAINTVIOLATIONS_HPP_)
#error "This file should only be included through ConstraintViolations.hpp"
#endif

#ifndef CONSTRAINTVIOLATIONS_IMPL_HPP_
#define CONSTRAINTVIOLATIONS_IMPL_HPP_

void ConstraintViolations::execute()
{
    compute_norms();
    write_norms();
}

std::pair<double, double> ConstraintViolations::get_norms()
{
    return std::make_pair(m_Ham_norm, m_Mom_norm);
}

void ConstraintViolations::compute_norms()
{
    m_Ham_norm = m_gr_amr->compute_norm(Interval(m_Ham_comp, m_Ham_comp),
                                      m_norm_exponent, m_dx_coarse);
    m_Mom_norm = m_gr_amr->compute_norm(m_Mom_comps, m_norm_exponent,m_dx_coarse);
}

void ConstraintViolations::write_norms() const
{
    int rank;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
    rank = 0;
#endif
    // only rank 0 does the write out
    if (rank == 0)
    {
        std::ofstream outfile;
        // overwrite file if this is the first timestep, otherwise append.
        if (m_time == m_dt)
        {
            outfile.open(m_filename);
        }
        else
        {
            outfile.open(m_filename, std::ios_base::app);
        }
        if (!outfile)
        {
            MayDay::Error(
                "ConstraintViolations::write_norms: error opening output file");
        }

        //Header data at first timestep
        if (m_time == m_dt)
        {
            outfile << '#' << std::setw(9) << "time";
            outfile << std::setw(20) << "L^" << m_norm_exponent << "_Ham";
            outfile << std::setw(20) << "L^" << m_norm_exponent << "_Mom\n";
        }

        // Now the data
        outfile << std::fixed << std::setw(10) << m_time;
        outfile << std::scientific << std::setprecision(10);
        outfile << std::setw(20) << m_Ham_norm;
        outfile << std::setw(20) << m_Mom_norm;
        outfile << std::endl;

        outfile.close();
    }
}

#endif /* CONSTRAINTVIOLATIONS_IMPL_HPP_ */
