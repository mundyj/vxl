#ifdef __GNUC__
#pragma implementation
#endif

//:
// \file

#include <vcl_sstream.h>
#include "vnl_io_bignum.h"
#include <vsl/vsl_binary_io.h>
#include <vsl/vsl_indent.h>
#include <vsl/vsl_binary_explicit_io.h>

//=================================================================================
//: Binary save self to stream.
void vsl_b_write(vsl_b_ostream & os, const vnl_bignum & p)
{
  const short io_version_no = 1;
  vsl_b_write(os, io_version_no);
  vcl_stringstream ss;
  ss << p;
  vsl_b_write(os, ss.str());
}

//=================================================================================
//: Binary load self from stream.
void vsl_b_read(vsl_b_istream &is, vnl_bignum & p)
{
  if (!is) return;
  short ver;
  vcl_string s;
  vcl_stringstream ss;
  vsl_b_read(is, ver);
  switch(ver)
  {
  case 1:
    vsl_b_read(is, s);
    ss.str(s);
    ss >> p;
    break;

  default:
    vcl_cerr << "I/O ERROR: vsl_b_read(vsl_b_istream&, vnl_bignum&) \n";
    vcl_cerr << "           Unknown version number "<< ver << "\n";
    is.is().clear(vcl_ios::badbit); // Set an unrecoverable IO error on stream
    return;
  }
}

//====================================================================================
//: Output a human readable summary to the stream
void vsl_print_summary(vcl_ostream & os,const vnl_bignum & p)
{
  os<<p;
}
