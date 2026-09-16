// Project: kepler-obake
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef KEPLER_SERIES_HPP
#define KEPLER_SERIES_HPP

#include <cstdint>

#include <obake/symbols.hpp>

#include "types.hpp"

namespace kepler
{

/// A Poisson series together with truncation metadata.
///
/// Wraps an obake power series and stores the total degree at which it
/// has been truncated, together with the set of *variables* (as opposed
/// to *parameters*) on which the truncation operates. Concrete Keplerian
/// series (z1, z2, ...) are constructed through the static factory
/// members, which encode the analytic formula of each function and
/// truncate it to a prescribed degree.
class series
{
public:
    /// Build the Poisson series of \f$z_1\f$ truncated at \p truncation_degree
    /// in the variables \f$X\f$ and \f$\bar{X}\f$.
    ///
    /// \f[
    ///   z_1 = \frac{i}{2}\,\bigl[X\Lambda^{-1} - \bar{X}\Lambda\bigr]\,
    ///         \sqrt{1 - \tfrac{1}{4} X\bar{X}}
    /// \f]
    [[nodiscard]] static series z1(::std::int64_t truncation_degree);

    /// Build the Poisson series of \f$z_2\f$ truncated at \p truncation_degree
    /// in the variables \f$X\f$ and \f$\bar{X}\f$.
    [[nodiscard]] static series z2(::std::int64_t truncation_degree);

    /// Const access to the underlying obake series.
    [[nodiscard]] const pser_t &get_series() const;

    /// The total degree at which the series was truncated in the variables.
    [[nodiscard]] ::std::int64_t get_truncation_degree() const;

    /// The set of variables on which truncation operates (e.g. {"X","Xc"}).
    [[nodiscard]] const obake::symbol_set &get_variables() const;

private:
    series(pser_t s, ::std::int64_t truncation_degree, obake::symbol_set variables);

    pser_t m_series;
    ::std::int64_t m_truncation_degree;
    obake::symbol_set m_variables;
};

namespace detail
{

/// Expand \f$\sqrt{1 - t}\f$ as a binomial series in
/// \f$t\f$ up to order \p truncation_degree.
[[nodiscard]] pser_t sqrt_one_minus_t(pser_t t, ::std::int64_t truncation_degree);
/// Expand \f$\sin t\f$ as a series in
/// \f$t\f$ up to order \p truncation_degree.
[[nodiscard]] pser_t sin(pser_t t, ::std::int64_t truncation_degree);
/// Expand \f$\cos t\f$ as a series in
/// \f$t\f$ up to order \p truncation_degree.
[[nodiscard]] pser_t cos(pser_t t, ::std::int64_t truncation_degree);

} // namespace detail

} // namespace kepler

#endif
