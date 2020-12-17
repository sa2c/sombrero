#include "communications.h"
#include "dirac.h"
#include "global.h"
#include "linear_algebra.h"
#include "random.h"
#include "suN.h"
#include "update.h"
#include <math.h>

void flat_source(spinor_field *s) {
  _ONE_SPINOR_FOR(s) { _spinor_one_f(*_SPINOR_PTR(s)); }
  // if(CID == 0){
  // 	_FIELD_AT(s,0)->c[0].c[0].re = 1.0;
  //}
  apply_BCs_on_spinor_field(s);

  start_sf_sendrecv(s);
  complete_sf_sendrecv(s);
}
