//------------------------------------------------------------------------------
// FileName    : bcal\bcam\camera_graph.h
// Author      : Kongbin Kang (kk@lems.brown.edu)
// Company     : Brown University
// Purpose     : a wapper class for provide some functionals
//               which are not supported by vcsl_graph.
// Date Of Creation: 3/25/2003
// Modification History :
// Date             Modifications
// 4/1/03           no longer use vcsl_graph. it use adjacent list to store edges
//------------------------------------------------------------------------------

#ifndef AFX_CAMERAGRAPH_H__8810B4C3_E6C7_42CE_8C7F_D9F8F201F6F6__INCLUDED_
#define AFX_CAMERAGRAPH_H__8810B4C3_E6C7_42CE_8C7F_D9F8F201F6F6__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <vcl_cassert.h>
#include <vcl_iostream.h>
#include <vcl_cstdlib.h>
#include <vcl_vector.h>
#include <vcsl/vcsl_spatial.h>
#include <vcsl/vcsl_spatial_sptr.h>
#include <vcsl/vcsl_spatial_transformation_sptr.h>


template<class ObjectNode, class CameraNode, class Trans>
class camera_graph
{
    protected:
      struct vertex_node{
      public:
        int id_;
        CameraNode *v_; // distinguish source with other node
        ObjectNode *s_;
      public:
        vertex_node(int id) 
        { 
          id_ = id;
          v_ = 0;
          s_ = 0;
        }
        
        ~vertex_node()
        {
          if(v_)
            delete v_;
          v_ = 0;
          
          if(s_)
            delete s_;
          s_ = 0;
        }
      };
      
      
      struct edge_node{ 
        vertex_node* v_; // recode which vertex attached with this edge
        Trans* e_;
      public:
        edge_node(vertex_node* v) { v_ = v; e_ = new Trans;}
        virtual ~edge_node() {  if(e_) delete e_;   }
      };
      
        
      protected:
        int init_graph()
        {
          source_ = 0;
          num_vertice_ = 0;
          // add a source
          add_vertex(0);
          
          source_ = vertice_[0]->s_;
          return 0;
        }
        
        vertex_node* malloc_new_node(bool is_source)
        {
          vertex_node *v = new vertex_node(num_vertice_++);
          
          if(is_source){
            v->s_ = new ObjectNode;
            v->v_ = 0;
          }
          else{
            v->s_ = 0;
            v->v_ = new CameraNode;
          }
          
          return v;
        }
        
  public:
    class iterator{
    public:
      iterator() :pos_ (0), _Ptr(0) {}
      iterator(vcl_vector<vertex_node*> *_P, int pos = 0) : _Ptr(_P), pos_(pos) {}
      iterator(const iterator& _X) : _Ptr(_X._Ptr), pos_(_X.pos_) {}
      CameraNode& operator*() const {return *((*_Ptr)[pos_]->v_); }
      CameraNode* operator->() const {return &(* *this); }

      iterator& operator=(iterator& _X)
      {
        pos_ = _X.pos_;
        _Ptr = _X._Ptr;
        return *this;
      }

      iterator& operator++() 
      {
        _Inc();
        return (*this); 
      }
      
      iterator operator++(int)
      {
        iterator _Tmp = *this;
        ++*this;
        return (_Tmp); 
      }
      
      iterator& operator--()
      {
        _Dec();
        return (*this); 
      }
      
      iterator operator--(int t)
      {
        const_iterator _Tmp = *this;
        --*this;
        return (_Tmp); 
      }
      
      bool operator==(const iterator& _X) const
      {return (_Ptr == _X._Ptr && pos_ = _X.pos_); }
      
      bool operator!=(const iterator& _X) const
      {return (!(*this == _X)); }
      
      void _Dec() { if(pos_>=0)  pos_--; }
      void _Inc() { pos_ ++; }
      
      vertex_node* node() const
      {return (*_Ptr)[pos]; }
      
      int node_id() const {return (*_ptr)[pos]->id_;}
      
    protected:
      vcl_vector<vertex_node*>* _Ptr;
      int pos_;
    };
    
    public: // constructor and de-constructor
      
      camera_graph()
      {
        init_graph();
      }
      
      virtual ~camera_graph()
      {
        erase_graph();
      }
      
    public: // operations
      ObjectNode* get_source() {  return source_;}
      int get_source_id() { return 0;}
      
      // add a new vertex and edge from the neighbour v to it
      // it return the position of the new vertex in the array
      int add_vertex(int neighbour = 0)
      {
        // currently only edge from source can be added
        assert(neighbour == 0);
        
        // the neighbour should exist already
        assert(neighbour <= num_vertice_);
        
        vertex_node *v;
        if(source_){
          v = malloc_new_node(false); // allocate an ordinary vertex
        }
        else
          v = malloc_new_node(true); // allocate a source vertex
        
        vertice_.push_back(v);
        
        // add a neighbour list into edges
        vcl_vector<edge_node*> *neighbour_list = new vcl_vector<edge_node*>;
        edges_.push_back(neighbour_list); 
        
        // update old adjacent neighbour list
        if(source_){ // if not source vertex is added
          edge_node *e = new edge_node(v);
          edges_[neighbour]->push_back(e);
        }
        
        return v->id_;
      }
      
      // id is the same of it is position in array
      inline CameraNode* get_vertex(int i)
      {
        assert(i>0); 
        return vertice_[i]->v_;
      }

      // return a iterator pointer to first vertex
      iterator first()
      { iterator iter(&vertice_, 1); return iter; }
    
      // get edge from v1 to v2
      inline Trans* get_edge(int v1, int v2)
      {
        assert(v1 == 0 && v2 <= num_vertice_); // only from souce to camera is avaible
        vcl_vector<edge_node*>* plist = edges_[v1];
        for(int i=0; i < plist->size(); i++){
          edge_node *e = (*plist)[i];
          vertex_node *v = e->v_;
          assert(v != 0); // no single edge exist
          
          if(v->id_ == v2){
            return e->e_;
          }
          else
            continue;
        }
        
        return 0; // cannot find edge
      }
            
      int num_vertice() { return vertice_.size() - 1; }
      
      
      // for debug
      void print()
      {
        vcl_cerr<<"print graph\n";
        for(int i=0; i<num_vertice_; i++){
          vcl_cerr<<"vertex id is: "<<vertice_[i]->id_<<"\t v is: "<< \
            vertice_[i]->v_ <<"\t s is: "<<vertice_[i]->s_<<"\n";
        }
      }
      
     
        int erase_graph()
        {
          num_vertice_ = vertice_.size();
          assert(num_vertice_ == edges_.size());
          
          if(num_vertice == 0) // empty graph
            return 0; // no error
          
          for(int i=0; i<num_vertice_; i++){
            // delete each vertex
            if(vertice_[i])
              delete vertice_[i];
            vertice_[i] = 0;
            
            // delete each edge node
            if(edges_[i]){ 
              int list_length = edges_[i]->size();
              for(int j=0; j<list_length; j++)
                delete (*(edges_[i]))[j];
              
              delete edges_[i];
            }
            edges_[i] = 0;
          }
          vertice_.clear();
          edges_.clear();
          
          source_ = 0;
          
          return 0; // no error
        }

      private:
        vcl_vector<vertex_node*> vertice_;
        vcl_vector<vcl_vector<edge_node*>* > edges_; // adjacent neigbhour list
        ObjectNode* source_;
        int num_vertice_;

};

#endif // AFX_CAMERAGRAPH_H__8810B4C3_E6C7_42CE_8C7F_D9F8F201F6F6__INCLUDED_
