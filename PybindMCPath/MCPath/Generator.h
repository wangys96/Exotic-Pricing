/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

#include <ql/qldefines.hpp>
#ifdef BOOST_MSVC
#  include <ql/auto_link.hpp>
#endif
#include <iostream>
#include <pybind11.h>
#include <numpy.h>
#include <ql/methods/montecarlo/all.hpp>
#include <ql/termstructures/volatility/all.hpp>
#include <ql/termstructures/yield/all.hpp>
#include <ql/processes/all.hpp>
#include <ql/time/all.hpp>

#include <MyPathGenerator.h>

namespace py = pybind11;
using namespace QuantLib;

typedef enum IRType {
    FlatRateCurve, SpotRateCurve, ForwardRateCurve, DiscountFactorCurve
};
typedef enum VolType{
    FlatVolCurve, SpotVolCurve
};
typedef enum ProcType{
    BS, BSM
};
typedef enum DCType {
    A365, AA, A360, F360
};
typedef enum BarType {
    NoBarrier,ConstBarrier,NonConstBarrier
};

Date _ParseDate(py::tuple date) 
{
    return(Date((Day)date[0].cast<int>(),
                (Month)date[1].cast<int>(),
                (Year)date[2].cast<int>())
        );
}

DayCounter _MakeDC(int dc_enum)
{
    switch (dc_enum) {
    case A365:
        return(Actual365Fixed());
    case AA:
        return(ActualActual());
    case A360:
        return(Actual360());
    case F360:
        return(Thirty360());
    default:
        throw std::invalid_argument("DayCounter type is not surppoted.");
    }
}

std::vector<Date> _Date2Vec(Date today, py::array_t<int>& input) 
{
    auto r = input.unchecked<1>();
    std::vector<Date> result;
    for (int i = 0; i < r.shape(0); i++)
        result.push_back(today + r(i) * Days);
    return(result);
}

template<class OutType>
std::vector<OutType> _Data2Vec(py::array_t<double>& input) 
{
    auto r = input.unchecked<1>();
    std::vector<OutType> result;
    for (int i = 0; i < r.shape(0); i++)
        result.push_back((OutType)r(i));
    return(result);
}

Handle<YieldTermStructure> _MakeIRCurve(Date today, int ir_type, py::array_t<int>& ir_term, py::array_t<double>& ir_data, int ir_dc)
{
    std::vector<Date> ir_term_vec(_Date2Vec(today, ir_term));
    std::vector<double> ir_data_vec(_Data2Vec<double>(ir_data));
    DayCounter dc = _MakeDC(ir_dc);

    Handle<YieldTermStructure> ir_curve;
    auto ir_arr = ir_data.unchecked<1>();
    switch (ir_type) {
    case FlatRateCurve:
        //std::cout << "Making FlatForward "<< (Rate)ir_arr(0) <<" at "<< today << std::endl;
        ir_curve = Handle<YieldTermStructure>(
            ext::shared_ptr<YieldTermStructure>(
                new FlatForward(today, (Rate)ir_arr(0), dc)
                )
            );
        break;
    case SpotRateCurve:
        //std::cout << "Making SpotCurve " << (Rate)ir_arr(0) << " at " << today << std::endl;
        ir_curve = Handle<YieldTermStructure>(
            ext::shared_ptr<YieldTermStructure>(
                new InterpolatedZeroCurve<Linear>(ir_term_vec, ir_data_vec, dc, Linear())
                )
            );
        break;
    case ForwardRateCurve:
        //std::cout << "Making ForwardCurve " << (Rate)ir_arr(0) << " at " << today << std::endl;
        ir_curve = Handle<YieldTermStructure>(
            ext::shared_ptr<YieldTermStructure>(
                new InterpolatedForwardCurve<Linear>(ir_term_vec, ir_data_vec, dc, Linear())
                )
            );
        break;
    case DiscountFactorCurve:
        //std::cout << "Making DiscountCurve " << (Rate)ir_arr(0) << " at " << today << std::endl;
        ir_curve = Handle<YieldTermStructure>(
            ext::shared_ptr<YieldTermStructure>(
                new InterpolatedDiscountCurve<Linear>(ir_term_vec, ir_data_vec, dc, Linear())
                )
            );
        break;
    default:
        throw std::invalid_argument("IR curve type is not surppoted.");
    }
    return(ir_curve);
}

Handle<BlackVolTermStructure> _MakeVolCurve(Date today, int vol_type, py::array_t<int>& vol_term, py::array_t<double>& vol_data, int vol_dc) {
    std::vector<Date> vol_term_vec(_Date2Vec(today, vol_term));
    std::vector<double> vol_data_vec(_Data2Vec<double>(vol_data));
    DayCounter dc = _MakeDC(vol_dc);

    Handle<BlackVolTermStructure> vol_curve;
    auto vol_arr = vol_data.unchecked<1>();
    switch (vol_type) {
    case FlatVolCurve:
        //std::cout << "Making FlatVol " << (Volatility)vol_arr(0) << " at " << today << std::endl;
        vol_curve = Handle<BlackVolTermStructure>(
            ext::shared_ptr<BlackVolTermStructure>(
                new BlackConstantVol(today, NullCalendar(), (Volatility)vol_arr(0), dc)
                )
            );
        break;
    case SpotVolCurve:
        //std::cout << "Making VolCurve " << (Volatility)vol_arr(0) << " at " << today << std::endl;
        vol_curve = Handle<BlackVolTermStructure>(
            ext::shared_ptr<BlackVolTermStructure>(
                new BlackVarianceCurve(today,vol_term_vec, vol_data_vec, dc)
                )
            );
        break;
    default:
        throw std::invalid_argument("Vol curve type is not surppoted.");
    }
    return(vol_curve);
}

py::array_t<double> GeneratePath(py::tuple today, int num, int steps, double tenor,
        int ir_type,  py::array_t<int> ir_term,  py::array_t<double> ir_data,  int ir_dc,
        int d_type,   py::array_t<int> d_term,   py::array_t<double> d_data,   int d_dc,
        int vol_type, py::array_t<int> vol_term, py::array_t<double> vol_data, int vol_dc,
        int upout_type,   py::array_t<bool> upout_ob,   py::array_t<double> upout_barrier,
        int downout_type, py::array_t<bool> downout_ob, py::array_t<double> downout_barrier,
        int proc_type,  py::array_t<double> output_matrix,
        bool bb = true, int skip = 0, int seed = 42)      
{
    Date todayDate(_ParseDate(today));
    
    Handle<Quote> S0(ext::shared_ptr<Quote>(new SimpleQuote(1.0)));

    //std::cout << "Making IR Curve "<< std::endl;
    Handle<YieldTermStructure> ir_curve(_MakeIRCurve(todayDate,ir_type,ir_term,ir_data,ir_dc));
    //std::cout << "Making Vol Curve " << std::endl;
    Handle<BlackVolTermStructure> vol_curve(_MakeVolCurve(todayDate, vol_type, vol_term, vol_data, vol_dc));
    
    //std::cout << "Making Process " << std::endl;
    ext::shared_ptr<GeneralizedBlackScholesProcess> process;
    if (proc_type == BSM)
    {
        Handle<YieldTermStructure> d_curve = _MakeIRCurve(todayDate, d_type, d_term, d_data, d_dc);
        process = ext::shared_ptr<GeneralizedBlackScholesProcess>(
            new BlackScholesMertonProcess(S0, d_curve, ir_curve, vol_curve)
            );
    }
    else if (proc_type == BS)
    {
        process = ext::shared_ptr<GeneralizedBlackScholesProcess>(
            new BlackScholesProcess(S0, ir_curve, vol_curve)
            );
    }
    else
        throw std::invalid_argument("Process type is not surppoted.");

    typedef LowDiscrepancy::rsg_type RSGType;
    //typedef MyPathGenerator<LowDiscrepancy>::sample_type PathType;
    RSGType rsg = LowDiscrepancy::make_sequence_generator((Size)steps, seed);

    //std::cout << "Making Generator " << std::endl;
    MyPathGenerator<RSGType> generator(process, (Time)tenor, (Size)steps, rsg, bb);
    
    for (int i = 0; i < skip; i++)
        generator.gen_bm();
        //generator.next();

    #define CHECK_INTERRUPT(row)                                                                \
        if (row % 10000 == 0 && PyErr_CheckSignals() != 0)                                      \
        {                                                                                       \
            std::cout << "MC simulation stopped. " << row << " paths are made." << std::endl;   \
            break;                                                                              \
        } 

    #define COPY_PATH(METHOD_NAME,...)                                                          \
        for (ssize_t row = 0; row < num; row++)                                                 \
        {                                                                                       \
            CHECK_INTERRUPT(row)                                                                \
            generator.gen_bm();                                                                 \
            generator.METHOD_NAME(arr, row, ##__VA_ARGS__);                                     \
        }
    
    auto arr = output_matrix.mutable_unchecked<2>();
    auto arr_upout_ob = upout_ob.mutable_unchecked<1>();
    auto arr_upout_barrier = upout_barrier.mutable_unchecked<1>();
    auto arr_downout_ob = downout_ob.mutable_unchecked<1>();
    auto arr_downout_barrier = downout_barrier.mutable_unchecked<1>();
    double upout_b, downout_b;

    if (upout_type == ConstBarrier)
        upout_b = arr_upout_barrier(0);
    if (downout_type == ConstBarrier)
        downout_b = arr_downout_barrier(0);

    //No Early Stop
    if ((upout_type == NoBarrier) && (downout_type == NoBarrier))
        COPY_PATH(copy_next)

    //Down Out Stop Barrier
    else if ((upout_type == NoBarrier) && (downout_type == ConstBarrier))
        COPY_PATH(copy_next_downout, arr_upout_ob, upout_b, arr_downout_ob, downout_b)
    else if ((upout_type == NoBarrier) && (downout_type == NonConstBarrier))
        COPY_PATH(copy_next_downout, arr_upout_ob, arr_upout_barrier, arr_downout_ob, arr_downout_barrier)

    //Up Out Stop Barrier
    else if ((upout_type == ConstBarrier) && (downout_type == NoBarrier))
        COPY_PATH(copy_next_upout, arr_upout_ob, upout_b, arr_downout_ob, downout_b)
    else if ((upout_type == NonConstBarrier) && (downout_type == NoBarrier))
        COPY_PATH(copy_next_upout, arr_upout_ob, arr_upout_barrier, arr_downout_ob, arr_downout_barrier)
    
    //Dual Stop Barriers
    else if ((upout_type == ConstBarrier) && (downout_type == ConstBarrier))
        COPY_PATH(copy_next_dualout, arr_upout_ob, upout_b, arr_downout_ob, downout_b)
    else if ((upout_type == ConstBarrier) && (downout_type == NonConstBarrier))
        COPY_PATH(copy_next_dualout, arr_upout_ob, upout_b, arr_downout_ob, arr_downout_barrier)
    else if ((upout_type == NonConstBarrier) && (downout_type == ConstBarrier))
        COPY_PATH(copy_next_dualout, arr_upout_ob, arr_upout_barrier, arr_downout_ob, downout_b)
    else if ((upout_type == NonConstBarrier) && (downout_type == NonConstBarrier))
        COPY_PATH(copy_next_dualout, arr_upout_ob, arr_upout_barrier, arr_downout_ob, arr_downout_barrier)

    return(output_matrix);
}


void GenerateRS(int num, int steps, double tenor, py::array_t<double> output_matrix, bool bb=true, int skip = 0, int seed=42)
{
    typedef LowDiscrepancy::rsg_type RSGType;
    RSGType rsg = LowDiscrepancy::make_sequence_generator((Size)steps, seed);
    MyRandomSequenceGenerator<RSGType> generator((Time)tenor, (Size)steps, rsg, bb);
    auto arr = output_matrix.mutable_unchecked<2>();

    for (int i = 0; i < skip; i++)
        generator.gen_bm();

    for (ssize_t row = 0; row < num; row++)
    {
        CHECK_INTERRUPT(row)
        generator.copy_bm(arr, row);
    }
}