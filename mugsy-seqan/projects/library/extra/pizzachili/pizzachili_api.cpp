// This file is auto-generated.
// DO NOT EDIT MANUALLY! CHANGES WILL BE OVERWRITTEN!

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

//#include "pizzachili_api.h"
#include <seqan/index/pizzachili_api.h>

#ifndef PCINDEX_CPP
extern "C" {
#endif
#include "interface.h"
#ifndef PCINDEX_CPP
}
#endif

// Needed by the FM index.
extern "C"
int init_ds_ssort(int adist, int bs_ratio);

namespace SEQAN_NAMESPACE_MAIN {

char* PCINDEX_NAME::error_index(int e) {
    return ::error_index(e);
}

int PCINDEX_NAME::build_index(
    impl::uchar_t* text, ulong length, char* build_options, impl::index_t* index
) {
    return ::build_index(text, length, build_options, index);
}

int PCINDEX_NAME::save_index(impl::index_t index, char* filename) {
    return ::save_index(index, filename);
}

int PCINDEX_NAME::load_index(char* filename, impl::index_t* index) {
    return ::load_index(filename, index);
}

int PCINDEX_NAME::free_index(impl::index_t index) {
    return ::free_index(index);
}

int PCINDEX_NAME::index_size(impl::index_t index, impl::ulong_t* size) {
    return ::index_size(index, size);
}

int PCINDEX_NAME::count(
    impl::index_t index, impl::uchar_t* pattern, impl::ulong_t length, impl::ulong_t* numocc
) {
    return ::count(index, pattern, length, numocc);
}

int PCINDEX_NAME::locate(
    impl::index_t index,
    impl::uchar_t* pattern,
    impl::ulong_t length,
    impl::ulong_t** occ,
    impl::ulong_t* numocc
) {
    return ::locate(index, pattern, length, occ, numocc);
}

int PCINDEX_NAME::get_length(impl::index_t index, impl::ulong_t* length) {
    return ::get_length(index, length);
}

int PCINDEX_NAME::extract(
    impl::index_t index,
    impl::ulong_t from,
    impl::ulong_t to,
    impl::uchar_t** snippet,
    impl::ulong_t* snippet_length
) {
    return ::extract(index, from, to, snippet, snippet_length);
}

int PCINDEX_NAME::display(
    impl::index_t index,
    impl::uchar_t* pattern,
    impl::ulong_t length,
    impl::ulong_t numc, 
    impl::ulong_t* numocc,
    impl::uchar_t** snippet_text,
    impl::ulong_t** snippet_length
) {
    return ::display(index, pattern, length, numc, numocc, snippet_text, snippet_length);
}


#ifdef PCINDEX_FM

int PCINDEX_NAME::init_ds_ssort(int adist, int bs_ratio) {
    return ::init_ds_ssort(adist, bs_ratio);
}

#endif

} // namespace SEQAN_NAMESPACE_MAIN
