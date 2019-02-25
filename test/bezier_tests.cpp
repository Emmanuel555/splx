#include "catch2/catch.hpp"
#include "bezier.h"

TEST_CASE( "bezier test" ) {


  splx::Bezier<double, 3U> bez(3);
  REQUIRE(bez.m_a==3);


}
