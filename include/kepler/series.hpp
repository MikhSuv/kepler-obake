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

    [[nodiscard]] static series z3(::std::int64_t truncation_degree);

    [[nodiscard]] static series z3_z1z2(const series &z1, const series &z2, const series &z3);
    /// Const access to the underlying obake series.
    [[nodiscard]] const pser_t &get_series() const;

    /// The total degree at which the series was truncated in the variables.
    [[nodiscard]] ::std::int64_t get_truncation_degree() const;

    /// The set of variables on which truncation operates (e.g. {"X","Xc"}).
    [[nodiscard]] const obake::symbol_set &get_variables() const;

    /// Save the series to the file \p path as a Boost text archive.
    ///
    /// Writes the object through a Boost serialization text archive.
    /// Throws a ::std::runtime_error if the file cannot be opened or written.
    void save(const ::std::string &path) const;

    /// Load a series from the file \p path as a Boost text archive.
    ///
    /// Reads the object back and returns it. Throws a ::std::runtime_error
    /// if the file cannot be opened or was written by an incompatible
    /// serialization version.
    [[nodiscard]] static series load(const ::std::string &path);

    /// Construct an empty series.
    ///
    /// Used for deserialization: the series is subsequently populated by
    /// reading from a Boost serialization archive.
    series() = default;

private:
    series(const pser_t &s, ::std::int64_t truncation_degree, const obake::symbol_set &variables);

    pser_t m_series;
    ::std::int64_t m_truncation_degree;
    obake::symbol_set m_variables;

    // Serialisation.
    friend class ::boost::serialization::access;

    /// The version of the serialization format used by this class.
    static inline constexpr unsigned s11n_version = 0u;

    /// Serialize the object into a Boost archive \p ar.
    ///
    /// Writes the format version, the truncation degree, the set of
    /// variables and the underlying obake power series (which, in turn,
    /// serializes its own symbol set, truncation policy and all its
    /// terms) into the archive.
    template <class Archive> void save(Archive &ar, unsigned) const
    {
        ar << s11n_version;
        ar << m_truncation_degree;
        ar << m_variables;
        ar << m_series;
    }

    /// Deserialize the object from a Boost archive \p ar.
    ///
    /// Throws if the archive was written with an incompatible version
    /// of the serialization format.
    template <class Archive> void load(Archive &ar, unsigned)
    {
        unsigned version;
        ar >> version;
        if (version != s11n_version) {
            throw ::std::runtime_error("Incompatible kepler::series serialization version" +
                                       ::std::to_string(version));
        }
        ar >> m_truncation_degree;
        ar >> m_variables;
        ar >> m_series;
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

namespace detail
{

/// Expand \f$\sqrt{1 - t}\f$ as a binomial series in
/// \f$t\f$ up to order \p truncation_degree.
[[nodiscard]] pser_t sqrt_one_minus_t(const pser_t &t, ::std::int64_t truncation_degree);
/// Expand \f$\sin t\f$ as a series in
/// \f$t\f$ up to order \p truncation_degree.
[[nodiscard]] pser_t sin(const pser_t &t, ::std::int64_t truncation_degree);
/// Expand \f$\cos t\f$ as a series in
/// \f$t\f$ up to order \p truncation_degree.
[[nodiscard]] pser_t cos(const pser_t &t, ::std::int64_t truncation_degree);

} // namespace detail

} // namespace kepler

#endif
