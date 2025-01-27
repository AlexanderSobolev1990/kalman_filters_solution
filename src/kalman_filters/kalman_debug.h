//------------------------------------------------------------------------------
///
/// \file       kalman_filter_debug.h
/// \brief      Отладочный ключ для фильтров Калмана
/// \date       25.04.21 - создан, 11.12.24 - review
/// \author     Соболев А.А.
/// \addtogroup kalman_filters
/// \{
///

#ifndef KALMAN_FILTER_DEBUG_H
#define KALMAN_FILTER_DEBUG_H

//------------------------------------------------------------------------------
#ifdef DEBUG_KALMAN
#undef DEBUG_KALMAN
#endif
// #define DEBUG_KALMAN // Включение отладочной печати матриц фильтра в консоль

#endif // KALMAN_FILTER_DEBUG_H
/// \}
