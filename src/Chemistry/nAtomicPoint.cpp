#include <chemistry.h>

N::ProteinBank::AtomicPoint:: AtomicPoint(void)
{
  Values [ 0 ] = 0.0 ;
  Values [ 1 ] = 0.0 ;
  Values [ 2 ] = 0.0 ;
}

N::ProteinBank::AtomicPoint:: AtomicPoint(float x,float y,float z)
{
  Values [ 0 ] = x ;
  Values [ 1 ] = y ;
  Values [ 2 ] = z ;
}

N::ProteinBank::AtomicPoint::~AtomicPoint (void)
{
}

float N::ProteinBank::AtomicPoint::x(void) const
{
  return Values [ 0 ] ;
}

float N::ProteinBank::AtomicPoint::y(void) const
{
  return Values [ 1 ] ;
}

float N::ProteinBank::AtomicPoint::z(void) const
{
  return Values [ 2 ] ;
}

const float * N::ProteinBank::AtomicPoint::values(void) const
{
  return Values ;
}

void N::ProteinBank::AtomicPoint::setX(float x)
{
  Values [ 0 ] = x ;
}

void N::ProteinBank::AtomicPoint::setY(float y)
{
  Values [ 1 ] = y ;
}

void N::ProteinBank::AtomicPoint::setZ(float z)
{
  Values [ 2 ] = z ;
}
