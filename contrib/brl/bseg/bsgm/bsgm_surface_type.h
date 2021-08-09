// This is//terra/bsgm_surface_type_image.h
#ifndef bsgm_surface_type_image_h
#define bsgm_surface_type_image_h

//:
// \file
// \brief A class to store surface_type probabilities
// \author J.L. Mundy
// \date August 4, 2021
//
// \verbatim
//  Modifications
//   <none yet>
// \endverbatim
#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <limits>
#include <string>
#include <vil/vil_image_view.h>
class bsgm_surface_type
{
 public:
  enum stype { NO_DATA, INVALID_DATA, SHADOW, SHADOW_STEP, NO_SURFACE_TYPE};
  enum source { RECTIFIED_TARGET, DSM, NO_SOURCE};

 bsgm_surface_type():ni_(0), nj_(0){init_type_names();}

 bsgm_surface_type(source s, size_t ni, size_t nj):source_(s), ni_(ni), nj_(nj){init_type_names(); init_type_images();}

  //: from list of surface_type images
 bsgm_surface_type(source const& s, std::map<stype, vil_image_view<float> > const& type_images): source_(s), type_images_(type_images){init_type_names();
    ni_ = type_images_[NO_DATA].ni();nj_ = type_images_[NO_DATA].nj();}

  //:load from tif files
 bsgm_surface_type(std::string const& directory) { this->load_surface_types(directory); }

 //: set type image layer
 bool set_type_image(stype type, vil_image_view<float> const& type_image){
   if((type_image.ni() != ni_) || (type_image.nj() != nj_))
     return false;
   type_images_[type] = type_image;
   return true;
 }
 //: get type probability (set as well)
 float& p(size_t i, size_t j, stype type){ return type_images_[type](i, j);}

 //: apply a bool image to set probabilites to 1.0f == true, 0.0f == false
 bool apply(vil_image_view<bool> const& mask, stype type);

 //: apply a probability image to set probabilites 
 bool apply(vil_image_view<float> const & prob, stype type);
 

 //: map string to surface_type index
 stype type_from_string(std::string const& type_string) {
    for(std::map<stype, std::string>::iterator cit = type_names_.begin();
        cit != type_names_.end(); ++cit){
      if(cit->second == type_string)
        return cit->first;
    }
    return NO_SURFACE_TYPE;
 }
  //: map surface_type index to string
  std::string type_to_string(stype const& type){return type_names_[type];}

  source source_from_string(std::string const& source_str) {
    if (source_str == "rectified_target") return RECTIFIED_TARGET;
    else if (source_str == "DSM") return DSM;
    return NO_SOURCE;
  }
  std::string source_to_string(source const& src) {
    if (src == RECTIFIED_TARGET)
      return "rectified_target";
    else if (src == DSM)
      return "DSM";
    return "no_source";
  }
  bool load_surface_types(std::string const& path);

  bool save_surface_types(std::string const& path);

  //: accessors
  size_t ni(){return ni_;}
  size_t nj(){return nj_;}
  source source_id() const {return source_;}
  std::vector<std::string> defined_types() const{
    std::vector<std::string> ret;
    for(std::map<stype, std::string>::const_iterator nit = type_names_.begin();
        nit != type_names_.end(); nit++) ret.push_back(nit->second);
    return ret;
  }
  vil_image_view<float>& type_image(std::string const& type_name){
    return type_images_[type_from_string(type_name)];
  }
  vil_image_view<float>& type_image(stype type){
    return type_images_[type];
  }
    private:
  // internal methods
  void init_type_names(){
    type_names_[NO_DATA] = "no_data";
    type_names_[INVALID_DATA] = "invalid_data";
    type_names_[SHADOW] = "shadow";
    type_names_[SHADOW_STEP] = "shadow_step";
  }
  void init_type_images(){
    for(std::map<stype, std::string>::iterator nit = type_names_.begin();
        nit != type_names_.end(); ++nit){
      type_images_[nit->first] = vil_image_view<float>(ni_, nj_);
      type_images_[nit->first].fill(0.0f);
    }
  }
  // members
  source source_; 
  size_t ni_;
  size_t nj_;
  std::map<stype, std::string> type_names_;
  std::map<stype, vil_image_view<float> > type_images_;
};
#endif
