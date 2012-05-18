/*
 * PairRange.hpp
 *
 *  Created on: May 16, 2012
 *      Author: david
 *
 * Provides foreach capability for an std::pair of iterators.
 *
 * http://stackoverflow.com/questions/6167598/why-was-pair-range-access-removed-from-c11
 *
 * Usage:
 *  std::multimap<int,int> mm;
 *  ...
 *  for (auto& p : as_range(mm.equal_range(42))) {
 *     ...
 *  }
 *
 */

#ifndef PAIRRANGE_HPP_
#define PAIRRANGE_HPP_

#include <boost/graph/adjacency_list.hpp>

//namespace boost_graph_helpers
//{

namespace detail
{
template<typename T, typename IndexMap = boost::identity_property_map>
    class ptr_vector_property_map
        : public boost::put_get_helper<
              typename std::iterator_traits<typename std::vector<T>::iterator>::reference,
              ptr_vector_property_map<T,IndexMap> >
    {
    public:
        typedef typename boost::property_traits<IndexMap>::key_type  key_type;
        typedef T value_type;
        typedef typename std::iterator_traits<typename std::vector<T>::iterator>::reference reference;
        typedef boost::lvalue_property_map_tag category;

        ptr_vector_property_map(std::vector<T>* data, const IndexMap& index = IndexMap())
        : store(data), index(index)
        {}

        typename std::vector<T>::iterator storage_begin() { return store->begin(); }

        typename std::vector<T>::iterator storage_end() { return store->end(); }

        typename std::vector<T>::const_iterator storage_begin() const { return store->begin();  }

        typename std::vector<T>::const_iterator storage_end() const { return store->end(); }

        IndexMap&       get_index_map()       { return index; }
        const IndexMap& get_index_map() const { return index; }

    public:
        // Copy ctor absent, default semantics is OK.
        // Assignment operator absent, default semantics is OK.
        // CONSIDER: not sure that assignment to 'index' is correct.

        reference operator[](const key_type& v) const {
            typename boost::property_traits<IndexMap>::value_type i = get(index, v);
            if (static_cast<unsigned>(i) >= store->size()) {
                store->resize(i + 1, T());
            }
            return (*store)[i];
        }
    private:
        // Conceptually, we have a vector of infinite size. For practical
        // purposes, we start with an empty vector and grow it as needed.
        // Note that we cannot store pointer to vector here -- we cannot
        // store pointer to data, because if copy of property map resizes
        // the vector, the pointer to data will be invalidated.
        // I wonder if class 'pmap_ref' is simply needed.
        std::vector<T>* store;
        IndexMap index;
    };
}

template<typename T>
detail::ptr_vector_property_map<T> make_ptr_vector_property_map(std::vector<T>* v)
{
	return detail::ptr_vector_property_map<T>(v);
}

//}

#endif
