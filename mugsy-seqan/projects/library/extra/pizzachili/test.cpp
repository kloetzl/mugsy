#include <iostream>

#include <seqan/index/pizzachili_api.h>

using namespace std;
using namespace seqan;

#define MAKE_TEST(which, build_options) \
    inline void test_##which() { \
        cout << #which " test" << endl; \
        impl::index_t index; \
        impl::error_t e = PizzaChiliApi##which::build_index( \
            (unsigned char*)("Hallo, Welt!"), 11, build_options, &index \
        ); \
        \
        cout << e << ": " << PizzaChiliApi##which::error_index(e) << endl; \
        \
        if (e != 0) { \
            cout << endl; \
            return; \
        } \
        \
        impl::uchar_t* text; \
        impl::ulong_t len; \
        e = PizzaChiliApi##which::extract(index, 0, 11, &text, &len); \
        \
        if (e != 0) { \
            cout << e << ": " << PizzaChiliApi##which::error_index(e) << \
                endl << endl; \
            return; \
        } \
        cout << text << " (" << len << ")" << endl; \
        free(text); \
        \
        e = PizzaChiliApi##which::save_index(index, "indexdata"); \
        \
        if (e != 0) { \
            cout << e << ": " << PizzaChiliApi##which::error_index(e) << \
                endl << endl; \
            return; \
        } \
        \
        e = PizzaChiliApi##which::free_index(index); \
        if (e != 0) { \
            cout << e << ": " << PizzaChiliApi##which::error_index(e) << \
                endl << endl; \
            return; \
        } \
        e = PizzaChiliApi##which::load_index("indexdata", &index); \
        if (e != 0) { \
            cout << e << ": " << PizzaChiliApi##which::error_index(e) << \
                endl << endl; \
            return; \
        } \
        cout << endl; \
    }

MAKE_TEST(AF, "")
MAKE_TEST(CCSA, "")
MAKE_TEST(FM, "-f 100 -a 2")
MAKE_TEST(LZ, "")
MAKE_TEST(RSA, "")
//MAKE_TEST(RLFM, "")
MAKE_TEST(SA, "copy_text")
MAKE_TEST(SADA, "")
//MAKE_TEST(SSA, "")

int main() {
    test_AF();
    test_CCSA();
    test_FM();
    test_RSA();
    //test_RLFM();
    test_SA();
    test_SADA();
    //test_SSA();
}
