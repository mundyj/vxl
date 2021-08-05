// This is//terra/bsgm_category_image.h
#ifndef bsgm_category_image_h
#define bsgm_category_image_h

//:
// \file
// \brief A class to store category probabilities
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
class bsgm_category_image
{
 public:
  enum category { NO_DATA, INVALID_DATA, SHADOW, ROOF_OVERHANG, NO_CATEGORY};
  enum source { RECTIFIED_TARGET, DSM, NO_SOURCE};

  bsgm_category_image(){init_cat_names();}

  //: from list of category images
 bsgm_category_image(source const& s, std::map<category, vil_image_view<float> > const& categories): source_(s), categories_(categories){init_cat_names();}

  //:load from tif files
 bsgm_category_image(std::string const& directory) { this->load_categories(directory); }

  //: map string to category index
  category cat_from_string(std::string const& cat_string) {
    for(std::map<category, std::string>::iterator cit = category_names_.begin();
        cit != category_names_.end(); ++cit){
      if(cit->second == cat_string)
        return cit->first;
    }
    return NO_CATEGORY;
  }
  //: map category index to string
  std::string cat_to_string(category const& cat){return category_names_[cat];}

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
  bool load_categories(std::string const& path);

  bool save_categories(std::string const& path);

  //: accessors
  source source_enum() const {return source_;}
  
    private:
  // internal methods
  void init_cat_names(){
    category_names_[NO_DATA] = "no_data";
    category_names_[INVALID_DATA] = "invalid_data";
    category_names_[SHADOW] = "shadow";
    category_names_[ROOF_OVERHANG] = "roof_overhang";
  }
  // members
  source source_; 
  std::map<category, std::string> category_names_;
  std::map<category, vil_image_view<float> > categories_;
};
#endif
