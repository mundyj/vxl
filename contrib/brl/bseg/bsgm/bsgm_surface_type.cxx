#include "bsgm_surface_type.h"
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vul/vul_file.h>
#include <vul/vul_file_iterator.h>
#include <fstream>
#include <stdexcept>

bool bsgm_surface_type::load_surface_types(std::string const& directory) {
  if (!vul_file::is_directory(directory)) {
    std::string message = "category directory not accessable" + directory;
    std::cout << message << std::endl;
    return false;
  }
  std::string glob = "directory/*.tif";
  for (vul_file_iterator fn = glob; fn; ++fn) {
    std::string path = fn();
    vil_image_view<float> cat_img = vil_load(path.c_str());
    if (cat_img.ni() == 0) {
      std::string message = "can't load type image from " + path;
      std::cout << message << std::endl;
      return false;
    }
    path = vul_file::strip_directory(path);
    std::string type_str = vul_file::strip_extension(path);
    stype t = this->type_from_string(type_str);
    type_images_[t] = cat_img;
  }
  std::ifstream istr(directory + "/" + "source.txt");
  if (!istr) {
    std::cout << "source.txt missing from directory" << std::endl;
    return false;
    std::string source_str;
    istr >> source_str;
    source_ = this->source_from_string(source_str);
  }
  ni_ = type_images_[NO_DATA].ni();
  nj_ = type_images_[NO_DATA].nj();
  return true;
}

bool bsgm_surface_type::save_surface_types(std::string const& directory) {
  if (!vul_file::is_directory(directory)) {
    std::string message = "category directory not accessable" + directory;
    std::cout << message << std::endl;
    return false;
  }
  for (std::map<stype, vil_image_view<float> >::iterator cit = type_images_.begin();
    cit != type_images_.end(); ++cit) {
    std::string path = directory + "/" + type_to_string(cit->first) + ".tif";
    if (!vil_save(cit->second, path.c_str())) {
      std::cout << "Can't save type image to " << path << std::endl;
      return false;
    }
  }
  std::string source_path = directory + "/source.txt";
  std::ofstream ostr(source_path.c_str());
  if (!ostr) {
    std::cout << "Can't write to " << source_path << std::endl;
    return false;
  }
  ostr << source_to_string(source_) << std::endl;
  ostr.close();
  return true;
}

bool bsgm_surface_type::apply(vil_image_view<bool> const& mask, stype type){
  if(type_images_.count(type) == 0){
    std::cout << "specified type " << type << " does not exist" << std::endl;
    return false;
  }
  vil_image_view<float>& type_img = type_images_[type];
  if(type_img.ni()!=ni_ || type_img.nj() != nj_){
    std::cout << "mismatch in type image dimensions" << std::endl;
    return false;
  }
  for(size_t j = 0; j<nj_; ++j)
    for(size_t i = 0; i<ni_; ++i){
      bool m = mask(i,j);
      if(m) type_img(i,j) = 1.0f;
      else type_img(i,j) = 0.0f;
    }
  return true;
}
bool bsgm_surface_type::apply(vil_image_view<float> const& prob, stype type)
{
  if (type_images_.count(type) == 0)
  {
    std::cout << "specified type " << type << " does not exist" << std::endl;
    return false;
  }
  vil_image_view<float> & type_img = type_images_[type];
  if (type_img.ni() != ni_ || type_img.nj() != nj_)
  {
    std::cout << "mismatch in type image dimensions" << std::endl;
    return false;
  }
  for (size_t j = 0; j < nj_; ++j)
    for (size_t i = 0; i < ni_; ++i)
    {
      float pr = prob(i, j);
      type_img(i, j) = pr;
    }
  return true;
}
