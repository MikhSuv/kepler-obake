// Project: kepler-obake
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef KEPLER_SERIES_HPP
#define KEPLER_SERIES_HPP

#include <cstdint>
#include <stdexcept>
#include <string>

#include <boost/serialization/access.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/split_member.hpp>

#include <obake/symbols.hpp>

#include "kepler/types.hpp"

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
///
/// The underlying key is a packed monomial in the variables
/// \f$X\f$, \f$\bar{X}\f$ and \f$\Lambda\f$. The variable
/// \f$\Lambda\f$ (mean longitude) and its powers are deliberately *not*
/// included in the truncation: they act as parameters that are
/// substituted numerically at evaluation time, whereas truncation is
/// imposed on the amplitude variables \f$X\f$ and \f$\bar{X}\f$ only.
class series
{
public:
    /// Build the Poisson series of \f$z_1\f$ truncated at @p truncation_degree
    /// in the variables \f$X\f$ and \f$\bar{X}\f$.
    ///
    /// In Keplerian elements \f$z_1 = e \sin M\f$, which in the complex
    /// variables reads
    /// \f[
    ///   z_1 = \frac{i}{2}\,\bigl[X\Lambda^{-1} - \bar{X}\Lambda\bigr]\,
    ///         \sqrt{1 - \tfrac{1}{4} X\bar{X}}
    /// \f]
    ///
    /// @param truncation_degree the total degree in \f$X, \bar{X}\f$
    ///                          retained in the result.
    ///
    /// @return the truncated series of \f$z_1\f$.
    [[nodiscard]] static series z1(::std::int64_t truncation_degree);

    /// Build the Poisson series of \f$z_2\f$ truncated at @p truncation_degree
    /// in the variables \f$X\f$ and \f$\bar{X}\f$.
    ///
    /// In Keplerian elements \f$z_2 = e \cos M\f$, which in the complex
    /// variables reads
    /// \f[
    ///   z_2 = \frac{1}{2}\,\bigl[X\Lambda^{-1} + \bar{X}\Lambda\bigr]\,
    ///         \sqrt{1 - \tfrac{1}{4} X\bar{X}}
    /// \f]
    ///
    /// @param truncation_degree the total degree in \f$X, \bar{X}\f$
    ///                          retained in the result.
    ///
    /// @return the truncated series of \f$z_2\f$.
    [[nodiscard]] static series z2(::std::int64_t truncation_degree);

    /// Build the Poisson series of \f$z_3\f$ truncated at @p truncation_degree
    /// in the variables \f$z_1\f$ and \f$z_2\f$.
    ///
    /// \f$z_3\f$ is the solution of the implicit equation
    /// \f[
    ///   z_3 = z_1 \cos z_3 + z_2 \sin z_3,
    /// \f]
    /// that is, the equation of the centre \f$z_3 = \nu - M = e \sin \nu\f$
    /// linking the true anomaly \f$\nu\f$ to the mean anomaly \f$M\f$.
    /// The trigonometric functions of \f$z_1\f$ and \f$z_2\f$ are *not*
    /// substituted: the result is expressed in the formal variables
    /// \f$z_1, z_2\f$ and is meant to be fed to @ref z3_z1z2 together
    /// with the series of @ref z1 and @ref z2.
    ///
    /// The equation is solved by successive approximation,
    /// \f$w_1 = z_1\f$, \f$w_n = z_1 \cos w_{n-1} + z_2 \sin w_{n-1}\f$.
    /// The iteration contracts because
    /// \f$\bigl|\partial_w (z_1 \cos w + z_2 \sin w)\bigr| \le e < 1\f$,
    /// hence each step raises the order of accuracy by one; @p truncation_degree
    /// steps therefore deliver a series valid to that degree.
    ///
    /// @param truncation_degree the total degree in \f$z_1, z_2\f$
    ///                          retained in the result.
    ///
    /// @return the truncated series of \f$z_3\f$.
    [[nodiscard]] static series z3(::std::int64_t truncation_degree);

    /// Build the Poisson series of \f$z_4 = \sin z_3\f$, truncated at
    /// @p truncation_degree in the variable \f$z_3\f$.
    ///
    /// The result is expressed in the formal variable \f$z_3\f$ and is meant
    /// to be substituted with the series produced by @ref z3_z1z2.
    ///
    /// @param truncation_degree the total degree in \f$z_3\f$ retained in
    ///                          the result.
    ///
    /// @return the truncated series of \f$z_4\f$.
    [[nodiscard]] static series z4(::std::int64_t truncation_degree);

    /// Build the Poisson series of \f$z_5 = \cos z_3\f$, truncated at
    /// @p truncation_degree in the variable \f$z_3\f$.
    ///
    /// The result is expressed in the formal variable \f$z_3\f$ and is meant
    /// to be substituted with the series produced by @ref z3_z1z2.
    ///
    /// @param truncation_degree the total degree in \f$z_3\f$ retained in
    ///                          the result.
    ///
    /// @return the truncated series of \f$z_5\f$.
    [[nodiscard]] static series z5(::std::int64_t truncation_degree);

    /// Close the @ref z3 template by substituting the series of \f$z_1\f$ and
    /// \f$z_2\f$ into it.
    ///
    /// The formal variables \f$z_1, z_2\f$ appearing in the template built by
    /// @ref z3 are replaced by the corresponding expansions in the amplitude
    /// variables \f$X, \bar{X}\f$. Because every monomial of the template
    /// carries the same total degree in \f$(z_1, z_2)\f$ as in
    /// \f$(X, \bar{X})\f$, truncating before and after the substitution are
    /// equivalent, and the truncation degree of the template can safely be
    /// inherited from the incoming series.
    ///
    /// @param z1 the series of \f$z_1\f$, supplying the truncation degree.
    /// @param z2 the series of \f$z_2\f$.
    /// @param z3 the template of \f$z_3\f$ produced by @ref z3.
    ///
    /// @return the series of \f$z_3\f$ as a function of \f$X, \bar{X}\f$.
    [[nodiscard]] static series z3_z1z2(const series &z1, const series &z2, const series &z3);
    [[nodiscard]] static series z4z5_z3(const series &z3, const series &z45);

    /// Const access to the underlying obake series.
    ///
    /// @return the stored power series.
    [[nodiscard]] const pser_t &get_series() const;

    /// The total degree at which the series was truncated in the variables.
    ///
    /// @return the total truncation degree.
    [[nodiscard]] ::std::int64_t get_truncation_degree() const;

    /// The set of variables on which truncation operates (e.g. {"X","Xc"}).
    ///
    /// @return the truncation variables.
    [[nodiscard]] const ::obake::symbol_set &get_variables() const;

    /// Save the series to the file @p path as a Boost text archive.
    ///
    /// Writes the object through a Boost serialization text archive.
    ///
    /// @param path the file to write to.
    ///
    /// @throws ::std::runtime_error if the file cannot be opened or written.
    void save(const ::std::string &path) const;

    /// Load a series from the file @p path as a Boost text archive.
    ///
    /// Reads the object back and returns it.
    ///
    /// @param path the file to read from.
    ///
    /// @return the deserialized series.
    ///
    /// @throws ::std::runtime_error if the file cannot be opened or was
    ///                              written by an incompatible serialization
    ///                              version.
    [[nodiscard]] static series load(const ::std::string &path);

    /// Construct an empty series.
    ///
    /// Used for deserialization: the series is subsequently populated by
    /// reading from a Boost serialization archive.
    series() = default;

private:
    // NOTE: the payload is taken by value rather than by const reference so
    // that call sites can hand over a computed series without paying for a
    // deep copy. Passing a const reference here would silently degrade every
    // std::move() at the call site into a copy.
    series(pser_t ser, ::std::int64_t truncation_degree, ::obake::symbol_set variables);

    pser_t m_series;
    ::std::int64_t m_truncation_degree{0};
    ::obake::symbol_set m_variables;

    // Serialisation.
    // NOTE: save()/load() below are deliberately named after the Boost
    // serialization hooks, and are overloads of - not alternatives to - the
    // file-based save()/load() declared above. BOOST_SERIALIZATION_SPLIT_MEMBER
    // dispatches on those exact names.
    friend class ::boost::serialization::access;

    /// The version of the serialization format used by this class.
    static constexpr unsigned serialization_version = 0U;

    /// Serialize the object into a Boost archive @p arch.
    ///
    /// Writes the format version, the truncation degree, the set of
    /// variables and the underlying obake power series (which, in turn,
    /// serializes its own symbol set, truncation policy and all its
    /// terms) into the archive.
    ///
    /// @param arch the archive to write to.
    /// @param version the format version, ignored on writing.
    template <class Archive>
    void save(Archive &arch, unsigned /* version */) const
    {
        arch << serialization_version;
        arch << m_truncation_degree;
        arch << m_variables;
        arch << m_series;
    }

    /// Deserialize the object from a Boost archive @p arch.
    ///
    /// Reads the format version, the truncation degree, the set of variables
    /// and the underlying obake power series out of the archive.
    ///
    /// @param arch the archive to read from.
    /// @param version the format version, ignored on reading.
    ///
    /// @throws ::std::runtime_error if the archive was written with an
    ///                              incompatible version of the serialization
    ///                              format.
    template <class Archive>
    void load(Archive &arch, unsigned /* version */)
    {
        unsigned version{0U};
        arch >> version;
        if (version != serialization_version) {
            throw ::std::runtime_error("Incompatible kepler::series serialization version: "
                                       + ::std::to_string(version));
        }
        arch >> m_truncation_degree;
        arch >> m_variables;
        arch >> m_series;
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

namespace detail
{

/// Expand \f$\sqrt{1 - t}\f$ as a binomial series in @p base up to order
/// @p truncation_degree.
///
/// The dummy argument @p base stands for the symbol that will later be
/// substituted, so that the expansion can be performed before the
/// truncating substitution is known.
///
/// @param base the expansion variable.
/// @param truncation_degree the highest power of @p base to retain.
///
/// @return the truncated binomial expansion.
[[nodiscard]] pser_t sqrt_one_minus(const pser_t &base, ::std::int64_t truncation_degree);

/// Expand \f$\sin(\texttt{base})\f$ as a Taylor series up to order
/// @p truncation_degree.
///
/// @param base the expansion variable.
/// @param truncation_degree the highest power of @p base to retain.
///
/// @return the truncated Taylor expansion.
[[nodiscard]] pser_t sin(const pser_t &base, ::std::int64_t truncation_degree);

/// Expand \f$\cos(\texttt{base})\f$ as a Taylor series up to order
/// @p truncation_degree.
///
/// @param base the expansion variable.
/// @param truncation_degree the highest power of @p base to retain.
///
/// @return the truncated Taylor expansion.
[[nodiscard]] pser_t cos(const pser_t &base, ::std::int64_t truncation_degree);

} // namespace detail

} // namespace kepler

#endif
