//----------------------------------------------------------------------------------------------------------------------
///
/// \file       timing.h
/// \brief      Класс подсчета времени выполнения кода
/// \date       25.09.20 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#ifndef SPML_TIMING_H
#define SPML_TIMING_H

// System includes:
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <cassert>
#include <ratio>

namespace SPML /// Специальная библиотека программных модулей (СБ ПМ)
{
namespace Timing /// Функции работы со временем
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Класс для замера времени выполнения кода
///
class CTimeKeeper
{
public:
    void StartTimer();      ///< Запуск измерения
    void EndTimer();        ///< Текущее измерение завершено
    double TimeCur();       ///< Интервал текущий в секундах
    double TimeSumm();      ///< Интервал (суммарный) в секундах
    double TimePerOp();     ///< Интервал на 1 операцию в секундах
    int IntervalsNumber();  ///< Число замеренных интервалов

    CTimeKeeper(); ///< Конструктор

private:
//    std::chrono::steady_clock::duration durationSum_; ///< Интервал (суммарный)
    std::chrono::high_resolution_clock::duration durationSum_; ///< Интервал (суммарный)
    int intervalsNumber_;   ///< Число интервалов вошедших в суммарный интервал
    double timeSum_;        ///< Интервал (суммарный) в секундах
    double timePerOp_;      ///< Интервал на 1 операцию в секундах

//    std::chrono::steady_clock::time_point startCur_;    ///< Начало измерения (текущего)
//    std::chrono::steady_clock::time_point endCur_;      ///< Конец измерения (текущего)
//    std::chrono::steady_clock::duration durationCur_;   ///< Интервал (текущий)
    std::chrono::high_resolution_clock::time_point startCur_;    ///< Начало измерения (текущего)
    std::chrono::high_resolution_clock::time_point endCur_;      ///< Конец измерения (текущего)
    std::chrono::high_resolution_clock::duration durationCur_;   ///< Интервал (текущий)
    double timeCur_;    ///< Интервал (текущий) в секундах
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Текущее время
/// \return Строка времени в формате ЧЧ:ММ:СС.МММ
///
std::string Time_in_DD_MM_YY_HH_MM_SS_MMM_TZ();

}
}
#endif // SPML_TIMING_H
/// \}
