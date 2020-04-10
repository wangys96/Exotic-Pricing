/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

#include <ql/methods/montecarlo/brownianbridge.hpp>
#include <ql/stochasticprocess.hpp>
#include <pybind11.h>
#include <numpy.h>

namespace py = pybind11;
using namespace QuantLib;

typedef py::detail::unchecked_mutable_reference<double, 2i64> array2d_double;
typedef py::detail::unchecked_mutable_reference<double, 1i64> array1d_double;
typedef py::detail::unchecked_mutable_reference<bool, 1i64> array1d_bool;

template <class GSG>
class MyPathGenerator {
public:
    typedef Sample<Path> sample_type;
    // constructors
    MyPathGenerator(const ext::shared_ptr<StochasticProcess>&,
        Time length,
        Size timeSteps,
        const GSG& generator,
        bool brownianBridge);
    MyPathGenerator(const ext::shared_ptr<StochasticProcess>&,
        const TimeGrid& timeGrid,
        const GSG& generator,
        bool brownianBridge);
    //! \name inspectors
    //@{
    const sample_type& next() const;
    void gen_bm() const;
    void copy_bm  (array2d_double& arr, ssize_t& row) const;
    void copy_next(array2d_double& arr, ssize_t& row) const;
    void copy_term(array1d_double& drift, array1d_double& stoch) const;

    void copy_next_upout  (array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, double& upout_barrier,         array1d_bool& downout_ob, double& downout_barrier)         const;
    void copy_next_upout  (array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, array1d_double& upout_barrier, array1d_bool& downout_ob, array1d_double& downout_barrier) const;

    void copy_next_downout(array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, double& upout_barrier, array1d_bool& downout_ob, double& downout_barrier)         const;
    void copy_next_downout(array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, array1d_double& upout_barrier, array1d_bool& downout_ob, array1d_double& downout_barrier) const;

    void copy_next_dualout(array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, double& upout_barrier,         array1d_bool& downout_ob, double& downout_barrier)         const;
    void copy_next_dualout(array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, array1d_double& upout_barrier, array1d_bool& downout_ob, double& downout_barrier)         const;
    void copy_next_dualout(array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, double& upout_barrier,         array1d_bool& downout_ob, array1d_double& downout_barrier) const;
    void copy_next_dualout(array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, array1d_double& upout_barrier, array1d_bool& downout_ob, array1d_double& downout_barrier) const;

    const sample_type& antithetic() const;
    Size size() const { return dimension_; }
    const TimeGrid& timeGrid() const { return timeGrid_; }
    //@}
private:
    const sample_type& next(bool antithetic) const;
    bool brownianBridge_;
    GSG generator_;
    Size dimension_;
    TimeGrid timeGrid_;
    ext::shared_ptr<StochasticProcess1D> process_;
    mutable sample_type next_;
    mutable std::vector<Real> temp_;
    BrownianBridge bb_;
};

template <class GSG>
MyPathGenerator<GSG>::MyPathGenerator(
    const ext::shared_ptr<StochasticProcess>& process,
    Time length,
    Size timeSteps,
    const GSG& generator,
    bool brownianBridge)
    : brownianBridge_(brownianBridge), generator_(generator),
    dimension_(generator_.dimension()), timeGrid_(length, timeSteps),
    process_(ext::dynamic_pointer_cast<StochasticProcess1D>(process)),
    next_(Path(timeGrid_), 1.0), temp_(dimension_), bb_(timeGrid_) {
    QL_REQUIRE(dimension_ == timeSteps,
        "sequence generator dimensionality (" << dimension_
        << ") != timeSteps (" << timeSteps << ")");
}

template <class GSG>
MyPathGenerator<GSG>::MyPathGenerator(
    const ext::shared_ptr<StochasticProcess>& process,
    const TimeGrid& timeGrid,
    const GSG& generator,
    bool brownianBridge)
    : brownianBridge_(brownianBridge), generator_(generator),
    dimension_(generator_.dimension()), timeGrid_(timeGrid),
    process_(ext::dynamic_pointer_cast<StochasticProcess1D>(process)),
    next_(Path(timeGrid_), 1.0), temp_(dimension_), bb_(timeGrid_) {
    QL_REQUIRE(dimension_ == timeGrid_.size() - 1,
        "sequence generator dimensionality (" << dimension_
        << ") != timeSteps (" << timeGrid_.size() - 1 << ")");
}

template <class GSG>
const typename MyPathGenerator<GSG>::sample_type&
    MyPathGenerator<GSG>::next() const {
    return next(false);
}

template <class GSG>
const typename MyPathGenerator<GSG>::sample_type&
    MyPathGenerator<GSG>::antithetic() const {
    return next(true);
}

template <class GSG>
const typename MyPathGenerator<GSG>::sample_type&
    MyPathGenerator<GSG>::next(bool antithetic) const {

    typedef typename GSG::sample_type sequence_type;
    const sequence_type& sequence_ =
        antithetic ? generator_.lastSequence()
        : generator_.nextSequence();

    if (brownianBridge_) {
        bb_.transform(sequence_.value.begin(),
            sequence_.value.end(),
            temp_.begin());
    }
    else {
        std::copy(sequence_.value.begin(),
            sequence_.value.end(),
            temp_.begin());
    }

    next_.weight = sequence_.weight;

    Path& path = next_.value;
    path.front() = process_->x0();

    for (Size i = 1; i < path.length(); i++) {
        Time t = timeGrid_[i - 1];
        Time dt = timeGrid_.dt(i - 1);
        path[i] = process_->evolve(t, path[i - 1], dt,
            antithetic ? -temp_[i - 1] :
            temp_[i - 1]);
    }

    return next_;
}

//===================
// Custom Copy Path
//===================

template <class GSG>
void MyPathGenerator<GSG>::gen_bm() const
{
    typedef typename GSG::sample_type sequence_type;
    const sequence_type& sequence_ = generator_.nextSequence();

    if (brownianBridge_) {
        bb_.transform(sequence_.value.begin(),
            sequence_.value.end(),
            temp_.begin());
    }
    else {
        std::copy(sequence_.value.begin(),
            sequence_.value.end(),
            temp_.begin());
    }
    next_.weight = sequence_.weight;
}

template <class GSG>
void MyPathGenerator<GSG>::copy_term(array1d_double& drift, array1d_double& stoch) const
{
    Path& path = next_.value;
    Real last = 1;
    for (Size i = 1; i < path.length(); i++) {
        Time t = timeGrid_[i - 1];
        Time dt = timeGrid_.dt(i - 1);
        Real drift_value = process_->expectation(t, last, dt);
        Real stoch_value = process_->stdDeviation(t, last, dt);
        drift(i) = drift_value;
        stoch(i) = stoch_value;         
        last = process_->apply(drift_value, stoch_value*temp_[i - 1]);
    }
}

template <class GSG>
void MyPathGenerator<GSG>::copy_bm(array2d_double& arr, ssize_t& row) const
{
    typedef typename GSG::sample_type sequence_type;
    const sequence_type& sequence_ = generator_.nextSequence();

    if (brownianBridge_) {
        bb_.transform(sequence_.value.begin(),
            sequence_.value.end(),
            temp_.begin());
    }
    else {
        std::copy(sequence_.value.begin(),
            sequence_.value.end(),
            temp_.begin());
    }
    for (Size i = 1; i < next_.value.length(); i++)
        arr(row, i) = temp_[i-1];

    //next_.weight = sequence_.weight;
}

template <class GSG>
void MyPathGenerator<GSG>::copy_next(array2d_double& arr, ssize_t& row) const
{
    typedef typename GSG::sample_type sequence_type;
    const sequence_type& sequence_ = generator_.nextSequence();

    if (brownianBridge_) {
        bb_.transform(sequence_.value.begin(),
            sequence_.value.end(),
            temp_.begin());
    }
    else {
        std::copy(sequence_.value.begin(),
            sequence_.value.end(),
            temp_.begin());
    }

    next_.weight = sequence_.weight;

    Path& path = next_.value;
    Real last = 1;
    arr(row, 0) = 1;
    for (Size i = 1; i < path.length(); i++) {
        Time t = timeGrid_[i - 1];
        Time dt = timeGrid_.dt(i - 1);
        Real new_value = process_->evolve(t, last, dt, temp_[i - 1]);
        arr(row, i) = new_value;
        last = new_value;
    }
}

template <class GSG>
void MyPathGenerator<GSG>::copy_next_upout(array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, double& upout_barrier, array1d_bool& downout_ob, double& downout_barrier) const
{
    Path& path = next_.value;
    Real last = 1;
    arr(row, 0) = 1;
    for (Size i = 1; i < path.length(); i++) {
        Time t = timeGrid_[i - 1];
        Time dt = timeGrid_.dt(i - 1);
        Real new_value = process_->evolve(t, last, dt, temp_[i - 1]);
        arr(row, i) = new_value;
        if (upout_ob(i) && new_value >= upout_barrier)
            break;
        last = new_value;
    }
}

template <class GSG>
void MyPathGenerator<GSG>::copy_next_upout(array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, array1d_double& upout_barrier, array1d_bool& downout_ob, array1d_double& downout_barrier) const
{
    Path& path = next_.value;
    Real last = 1;
    arr(row, 0) = 1;
    for (Size i = 1; i < path.length(); i++) {
        Time t = timeGrid_[i - 1];
        Time dt = timeGrid_.dt(i - 1);
        Real new_value = process_->evolve(t, last, dt, temp_[i - 1]);
        arr(row, i) = new_value;
        if (upout_ob(i) && new_value >= upout_barrier(i))
            break;
        last = new_value;
    }
}

template <class GSG>
void MyPathGenerator<GSG>::copy_next_downout(array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, double& upout_barrier, array1d_bool& downout_ob, double& downout_barrier) const
{
    Path& path = next_.value;
    Real last = 1;
    arr(row, 0) = 1;
    for (Size i = 1; i < path.length(); i++) {
        Time t = timeGrid_[i - 1];
        Time dt = timeGrid_.dt(i - 1);
        Real new_value = process_->evolve(t, last, dt, temp_[i - 1]);
        arr(row, i) = new_value;
        if (downout_ob(i) && new_value < downout_barrier)
            break;
        last = new_value;
    }
}

template <class GSG>
void MyPathGenerator<GSG>::copy_next_downout(array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, array1d_double& upout_barrier, array1d_bool& downout_ob, array1d_double& downout_barrier) const
{
    Path& path = next_.value;
    Real last = 1;
    arr(row, 0) = 1;
    for (Size i = 1; i < path.length(); i++) {
        Time t = timeGrid_[i - 1];
        Time dt = timeGrid_.dt(i - 1);
        Real new_value = process_->evolve(t, last, dt, temp_[i - 1]);
        arr(row, i) = new_value;
        if (downout_ob(i) && new_value < downout_barrier(i))
            break;
        last = new_value;
    }
}

template <class GSG>
void MyPathGenerator<GSG>::copy_next_dualout(array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, array1d_double& upout_barrier, array1d_bool& downout_ob, double& downout_barrier) const
{
    Path& path = next_.value;
    Real last = 1;
    arr(row, 0) = 1;
    for (Size i = 1; i < path.length(); i++) {
        Time t = timeGrid_[i - 1];
        Time dt = timeGrid_.dt(i - 1);
        Real new_value = process_->evolve(t, last, dt, temp_[i - 1]);
        arr(row, i) = new_value;
        if ((upout_ob(i) && new_value >= upout_barrier(i)) || (downout_ob(i) && new_value < downout_barrier))
            break;
        last = new_value;
    }
}

template <class GSG>
void MyPathGenerator<GSG>::copy_next_dualout(array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, double& upout_barrier, array1d_bool& downout_ob, array1d_double& downout_barrier) const
{
    Path& path = next_.value;
    Real last = 1;
    arr(row, 0) = 1;
    for (Size i = 1; i < path.length(); i++) {
        Time t = timeGrid_[i - 1];
        Time dt = timeGrid_.dt(i - 1);
        Real new_value = process_->evolve(t, last, dt, temp_[i - 1]);
        arr(row, i) = new_value;
        if ((upout_ob(i) && new_value >= upout_barrier) || (downout_ob(i) && new_value < downout_barrier(i)))
            break;
        last = new_value;
    }
}

template <class GSG>
void MyPathGenerator<GSG>::copy_next_dualout(array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, double& upout_barrier, array1d_bool& downout_ob, double& downout_barrier) const
{
    Path& path = next_.value;
    Real last = 1;
    arr(row, 0) = 1;
    for (Size i = 1; i < path.length(); i++) {
        Time t = timeGrid_[i - 1];
        Time dt = timeGrid_.dt(i - 1);
        Real new_value = process_->evolve(t, last, dt, temp_[i - 1]);
        arr(row, i) = new_value;
        if ((upout_ob(i) && new_value >= upout_barrier) || (downout_ob(i) && new_value < downout_barrier))
            break;
        last = new_value;
    }
}

template <class GSG>
void MyPathGenerator<GSG>::copy_next_dualout(array2d_double& arr, ssize_t& row, array1d_bool& upout_ob, array1d_double& upout_barrier, array1d_bool& downout_ob, array1d_double& downout_barrier) const
{
    Path& path = next_.value;
    Real last = 1;
    arr(row, 0) = 1;
    for (Size i = 1; i < path.length(); i++) {
        Time t = timeGrid_[i - 1];
        Time dt = timeGrid_.dt(i - 1);
        Real new_value = process_->evolve(t, last, dt, temp_[i - 1]);
        arr(row, i) = new_value;
        if ((upout_ob(i) && new_value >= upout_barrier(i)) || (downout_ob(i) && new_value < downout_barrier(i)))
            break;
        last = new_value;
    }
}


//===================
// Custom RSG
//===================
/*
template <class GSG>
class MyRandomSequenceGenerator {
public:
    typedef Sample<Path> sample_type;
    // constructors
    MyRandomSequenceGenerator(
        Time length,
        Size timeSteps,
        const GSG& generator,
        bool brownianBridge);
    const TimeGrid& timeGrid() const { return timeGrid_; }
    void copy_next(array2d_double& input_matrix, ssize_t& row) const;
    //@}
private:    
    bool brownianBridge_;
    GSG generator_;
    Size dimension_;
    TimeGrid timeGrid_;
    mutable sample_type next_;
    mutable std::vector<Real> temp_;
    BrownianBridge bb_;
};

template <class GSG>
MyRandomSequenceGenerator<GSG>::MyRandomSequenceGenerator( 
    Time length,
    Size steps,
    const GSG& generator,
    bool brownianBridge)
    : brownianBridge_(brownianBridge), generator_(generator),
    dimension_(generator_.dimension()), timeGrid_(length, steps),
    next_(Path(timeGrid_), 1.0), temp_(dimension_), bb_(timeGrid_) {
    QL_REQUIRE(dimension_ == steps,
        "sequence generator dimensionality (" << dimension_
        << ") != steps (" << steps << ")");
}

template <class GSG>
void MyRandomSequenceGenerator<GSG>::copy_next(array2d_double& input_matrix, ssize_t& row) const
{
    typedef typename GSG::sample_type sequence_type;
    const sequence_type& sequence_ = generator_.nextSequence();

    if (brownianBridge_) {
        bb_.transform(sequence_.value.begin(),
            sequence_.value.end(),
            temp_.begin());
    }
    else {
        std::copy(sequence_.value.begin(),
            sequence_.value.end(),
            temp_.begin());
    }
    for (Size i = 1; i < next_.value.length(); i++)
        input_matrix(row, i) = temp_[i - 1];
}
*/

template <class GSG>
class MyRandomSequenceGenerator {
public:
    typedef Sample<Path> sample_type;
    // constructors
    MyRandomSequenceGenerator(
        Time length,
        Size timeSteps,
        const GSG& generator,
        bool brownianBridge);

    void gen_bm() const;
    void copy_bm(array2d_double& arr, ssize_t& row) const;
    Size size() const { return dimension_; }
    const TimeGrid& timeGrid() const { return timeGrid_; }

private:
    bool brownianBridge_;
    GSG generator_;
    Size dimension_;
    TimeGrid timeGrid_;
    mutable sample_type next_;
    mutable std::vector<Real> temp_;
    BrownianBridge bb_;
};

template <class GSG>
MyRandomSequenceGenerator<GSG>::MyRandomSequenceGenerator(
    Time length,
    Size timeSteps,
    const GSG& generator,
    bool brownianBridge)
    : brownianBridge_(brownianBridge), generator_(generator),
    dimension_(generator_.dimension()), timeGrid_(length, timeSteps),
    next_(Path(timeGrid_), 1.0), temp_(dimension_), bb_(timeGrid_) {
    QL_REQUIRE(dimension_ == timeSteps,
        "sequence generator dimensionality (" << dimension_
        << ") != timeSteps (" << timeSteps << ")");
}

template <class GSG>
void MyRandomSequenceGenerator<GSG>::copy_bm(array2d_double& arr, ssize_t& row) const
{
    typedef typename GSG::sample_type sequence_type;
    const sequence_type& sequence_ = generator_.nextSequence();

    if (brownianBridge_) {
        bb_.transform(sequence_.value.begin(),
            sequence_.value.end(),
            temp_.begin());
    }
    else {
        std::copy(sequence_.value.begin(),
            sequence_.value.end(),
            temp_.begin());
    }
    for (Size i = 1; i < next_.value.length(); i++)
        arr(row, i) = temp_[i - 1];
}

template <class GSG>
void MyRandomSequenceGenerator<GSG>::gen_bm() const
{
    typedef typename GSG::sample_type sequence_type;
    const sequence_type& sequence_ = generator_.nextSequence();

    if (brownianBridge_) {
        bb_.transform(sequence_.value.begin(),
            sequence_.value.end(),
            temp_.begin());
    }
    else {
        std::copy(sequence_.value.begin(),
            sequence_.value.end(),
            temp_.begin());
    }
}