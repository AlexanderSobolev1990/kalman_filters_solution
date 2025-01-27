//----------------------------------------------------------------------------------------------------------------------
///
/// \file       radars.h
/// \brief      Структура определения позиции РЛС
/// \date       06.11.19 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#ifndef RADARS_H
#define RADARS_H

// SPML includes:
#include <geodesy.h>

namespace SPML /// Специальная библиотека программных модулей (СБПМ)
{
namespace Radar /// Позиции РЛС
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Информация об одном азимутальном секторе контроля изделия
/// \details Включает координаты, расстояние между приёмником и передатчиком, азимут от приёмника на передатчик,
/// центральный азимут сектора приёмника и т.д.
///
class CRadarPosition
{
public:    
    std::string Name() const; ///< Название сектора
    int Nsector() const; ///< Номер сектора обзора изделия
    double ACS() const; ///< Центральный абсолютный азимут сектора приемника
    double TransmitterLatitude() const; ///< Широта передатчика
    double TransmitterLongitude() const; ///< Долгота передатчика
    double RecieverLatitude() const; ///< Широта приемника
    double RecieverLongitude() const; ///< Долгота приемника
    double AzRT() const; ///< Азимут от приемника на передатчик
    double RRT() const; ///< Расстояние между приемником и передатчиком по Земле

    CRadarPosition(); ///< Конструктор по умолчанию

    ///
    /// \brief Параметрический конструктор
    /// \details Единицы измерений - согласно rangeUnit, angleUnit
    /// \param name - имя сектора
    /// \param nsector - номер сектора
    /// \param centralAzimuthOfSector - центральный азимут сектора
    /// \param transmitterLatitude - широта передатчика
    /// \param transmitterLongitude - долгота передатчика
    /// \param recieverLatitude - широта приемника
    /// \param recieverLongitude - долгота приемника
    /// \param ellipsoid - земной эллипсоид для расчета расстояния и азимута между позициями
    /// \param rangeUnit - единицы измерения расстояния
    /// \param angleUnit - единицы измерения углов (азимута)
    ///
    CRadarPosition( std::string name, int nsector, double centralAzimuthOfSector,
        double transmitterLatitude, double transmitterLongitude, double recieverLatitude, double recieverLongitude,
        const SPML::Geodesy::CEllipsoid &ellipsoid, const SPML::Units::TRangeUnit &rangeUnit, const SPML::Units::TAngleUnit &angleUnit );

private:
    std::string name_; ///< Название сектора
    int nsector_; ///< Номер сектора обзора изделия
    double acs_; ///< Центральный абсолютный азимут сектора приемника
    double transmitterLatitude_; ///< Широта передатчика
    double transmitterLongitude_; ///< Долгота передатчика
    double recieverLatitude_; ///< Широта приемника
    double recieverLongitude_; ///< Долгота приемника
    double azrt_; ///< Азимут от приемника на передатчик
    double rrt_; ///< Расстояние между приемником и передатчиком по Земле
};

} // end namespace Radar
} // end namespace SPML
#endif // RADARS_H
/// \}
