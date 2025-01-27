//----------------------------------------------------------------------------------------------------------------------
///
/// \file       radars.cpp
/// \brief      Все что связано с ?
/// \date       06.11.19 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#include <radars.h>

namespace SPML /// Специальная библиотека программных модулей (СБПМ)
{
namespace Radar /// Позиции РЛС
{

std::string CRadarPosition::Name() const { return name_; }

int CRadarPosition::Nsector() const { return nsector_; }

double CRadarPosition::ACS() const { return acs_; }

double CRadarPosition::TransmitterLatitude() const { return transmitterLatitude_; }

double CRadarPosition::TransmitterLongitude() const { return transmitterLongitude_; }

double CRadarPosition::RecieverLatitude() const { return recieverLatitude_; }

double CRadarPosition::RecieverLongitude() const { return recieverLongitude_; }

double CRadarPosition::AzRT() const { return azrt_; }

double CRadarPosition::RRT() const { return rrt_; }

CRadarPosition::CRadarPosition(){}

CRadarPosition::CRadarPosition( std::string name, int nsector, double centralAzimuthOfSector,
    double transmitterLatitude, double transmitterLongitude, double recieverLatitude, double recieverLongitude,
    const SPML::Geodesy::CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit )
{
    name_ = name;
    nsector_ = nsector;
    acs_ = centralAzimuthOfSector;
    transmitterLatitude_ = transmitterLatitude;
    transmitterLongitude_ = transmitterLongitude;
    recieverLatitude_ = recieverLatitude;
    recieverLongitude_ = recieverLongitude;
    SPML::Geodesy::GEOtoRAD( ellipsoid, rangeUnit, angleUnit, recieverLatitude, recieverLongitude,
        transmitterLatitude, transmitterLongitude, rrt_, azrt_ );
}

}
}
/// \}
