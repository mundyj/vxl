#include <vcl_iostream.h>
#include <vcl_fstream.h>
#include <vcl_utility.h>

#include <vnl/vnl_test.h>
#include <vnl/vnl_bignum.h>

#include <vnl/io/vnl_io_bignum.h>

void test_bignum_io()
{
  vcl_cout << "**************\n";
  vcl_cout << "test_bignum_io\n";
  vcl_cout << "**************\n";

  vnl_bignum nil(0), big(-3245444l), verybig("4235702934875938745092384750293845");

  vsl_b_ofstream bfs_out("vnl_bignum_test_io.bvl.tmp");
  TEST ("Created vnl_bignum_test_io.bvl.tmp for writing",
        (!bfs_out), false);
  vsl_b_write(bfs_out, nil);
  vsl_b_write(bfs_out, big);
  vsl_b_write(bfs_out, verybig);
  bfs_out.close();

  vnl_bignum r1, r2, r3;
  vsl_b_ifstream bfs_in("vnl_bignum_test_io.bvl.tmp");
  TEST ("Opened vnl_bignum_test_io.bvl.tmp for reading",
        (!bfs_in), false);
  vsl_b_read(bfs_in, r1);
  vsl_b_read(bfs_in, r2);
  vsl_b_read(bfs_in, r3);
  TEST ("Finished reading file successfully", (!bfs_in), false);
  bfs_in.close();

  TEST ("equality 0", nil, r1);
  TEST ("equality -3245444", big, r2);
  TEST ("equality 4235702934875938745092384750293845", verybig, r3);

  vcl_cout << "\n0 summary: ";
  vsl_print_summary(vcl_cout, nil);
  vcl_cout << "\n-3245444 summary: ";
  vsl_print_summary(vcl_cout, r2);
  vcl_cout << "\n4235702934875938745092384750293845 summary: ";
  vsl_print_summary(vcl_cout, verybig);
  vcl_cout << vcl_endl;
}

TESTMAIN(test_vector);
