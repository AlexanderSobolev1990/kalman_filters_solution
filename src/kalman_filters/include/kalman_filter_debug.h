//----------------------------------------------------------------------------------------------------------------------
///
/// \file       kalman_filter_debug.h
/// \brief      Отладочный ключ для фильтров Калмана
/// \date       25.04.21 - создан
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_DEBUG_H
#define KALMAN_FILTER_DEBUG_H

#ifdef DEBUG_KALMAN
#undef DEBUG_KALMAN
#endif
//#define DEBUG_KALMAN // Включение отладочной печати в консоль

#include <string>

namespace KalmanFilters /// Фильтры Калмана
{

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Возвращает строку, содержащую информацию о версии
/// \return Строка версии в формате DD-MM-YY-VV_COMMENTS, где DD - день, MM - месяц, YY - год, VV - версия, COMMENTS - комментарий(опционально)
///
std::string GetVersion();

}

#endif // KALMAN_FILTER_DEBUG_H
/// \}
