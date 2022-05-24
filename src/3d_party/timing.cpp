//----------------------------------------------------------------------------------------------------------------------
///
/// \file       timing.cpp
/// \brief      Класс подсчета времени выполнения кода
/// \date       25.09.20 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#include <timing.h>

namespace SPML /// Специальная библиотека программных модулей (СБ ПМ)
{
namespace Timing /// Функции работы со временем
{
//----------------------------------------------------------------------------------------------------------------------
CTimeKeeper::CTimeKeeper()
{
    intervalsNumber_ = 0;
    timeSum_ = 0;
    timePerOp_ = 0;
    timeCur_ = 0;
}

void CTimeKeeper::StartTimer()
{
    startCur_ = std::chrono::high_resolution_clock::now();
}

void CTimeKeeper::EndTimer()
{
    // Статистика по текущей операции
    endCur_ = std::chrono::high_resolution_clock::now();
    durationCur_ = endCur_ - startCur_;
    timeCur_ = double( durationCur_.count() ) * std::chrono::high_resolution_clock::period::num / std::chrono::high_resolution_clock::period::den;

    // Суммарная статистика
    intervalsNumber_++;
    timeSum_ += timeCur_;
    timePerOp_ = timeSum_ / static_cast<double>( intervalsNumber_ );
}

double CTimeKeeper::TimeCur()
{
    return timeCur_;
}

double CTimeKeeper::TimeSumm()
{
    return timeSum_;
}

double CTimeKeeper::TimePerOp()
{
    return timePerOp_;
}

int CTimeKeeper::IntervalsNumber()
{
    return intervalsNumber_;
}

//----------------------------------------------------------------------------------------------------------------------
std::string Time_in_DD_MM_YY_HH_MM_SS_MMM_TZ()
{    
    using namespace std::chrono;

    // get current time
    auto now = system_clock::now();

    // get number of milliseconds for the current second
    // (remainder after division into seconds)
    auto ms = duration_cast<milliseconds>( now.time_since_epoch() ) % 1000;

    // convert to std::time_t in order to convert to std::tm (broken-down time)
    auto timer = system_clock::to_time_t( now );

    // convert to broken-down time
    std::tm bdt = *std::localtime( &timer );

    std::ostringstream oss;

//    oss << std::put_time( &bt, "%d/%m/%y %H:%M:%S" );
    oss << std::put_time( &bdt, "%d-%m-%y %H:%M:%S" );

    auto gmtoffsec = bdt.tm_gmtoff; // Смещение отноистельно UTC в секундах
    std::string sign = "+";
    if( gmtoffsec < 0 ) {
        sign = "-";
    }
    auto gmtoffhour = std::trunc( gmtoffsec / 3600.0 );
    auto gmtoffmin = std::trunc( gmtoffsec / 60.0 ) - ( gmtoffhour * 60 );

    oss << '.' << std::setfill('0') << std::setw(3) << ms.count() << " UTC" << sign
        << std::setfill('0') << std::setw(2) << gmtoffhour
        << std::setfill('0') << std::setw(2) << gmtoffmin;

    return oss.str();
}

}
}
/// \}
