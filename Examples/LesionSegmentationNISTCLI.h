#ifndef __LesionSegmentationNISTCLI_h
#define __LesionSegmentationNISTCLI_h

#include "metaCommand.h"
#include "vtkImageData.h"
#include "vtkSmartPointer.h"
#include "itkFixedArray.h"
#include "itkLandmarkSpatialObject.h"
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

class LesionSegmentationNISTCLI : public MetaCommand
{

public:
  typedef itk::LandmarkSpatialObject< 3 >    SeedSpatialObjectType;
  typedef SeedSpatialObjectType::PointListType   PointListType;

  LesionSegmentationNISTCLI( int argc, char *argv[] ) : MetaCommand()
  {
    this->DisableDeprecatedWarnings();

    this->AddArgument("InputImage",false,"Input image to be segmented.");
    this->AddArgument("InputDICOMDir",false,"DICOM directory containing series of the Input image to be segmented.");
    this->AddArgument("OutputImage", false, "Output segmented image");
    this->AddArgument("OutputAIM", false, 
                      "Output segmented image in AIM format");
    this->AddArgument("OutputMesh", true, "Output segmented surface (STL filename expected)");
    this->AddArgument("OutputROI", false, "Write the ROI within which the segmentation will be confined to (for debugging purposes)");
    this->AddArgument("Visualize", false, "Visualize the input image and the segmented surface.", MetaCommand::BOOL, "0");
    this->AddArgument("Wireframe", false, "Visualize the input image and the segmented surface as a wireframe. Only valid if the Visualize flag is also enabled.", MetaCommand::BOOL, "0");
    this->AddArgument("IgnoreDirection", false, "Ignore the direction of the DICOM image", MetaCommand::BOOL, "0");
     this->AddArgument("PartSolid", false, "Default is to assume parameters for a solid lesion. Specify this if the lesion is part-solid.", MetaCommand::BOOL, "0");
    this->AddArgument("ROI", false,
                    "Bounds of the ROI if any, 6 parameters", MetaCommand::LIST);
    this->AddArgument("Sigma", false,
      "Manually specify sigma. This is an array with 3 values in physical units. This defaults to the maximum spacing in the dataset, if unspecified",
      MetaCommand::LIST);
    this->AddArgument("SeedFile", true,
                      "File with a list of seeds. At least one must be "
                      "specified.");
    this->AddArgument("MaximumRadius", false, "Maximum radius of the lesion in mm. This can be used as alternate way of specifying the bounds. You specify a seed and a value of say 20mm, if you know the lesion is smaller than 20mm..", MetaCommand::FLOAT, "30");
    
    
    if(!this->Parse(argc, argv))
      {
      // We can't invoke errors from constructors..
      exit(-1);
      }  
    }

  double *GetROI()
    {
    if (this->GetOptionWasSet("ROI"))
      {
      std::list< std::string > bounds = this->GetValueAsList("ROI");
      std::list< std::string >::const_iterator fit = bounds.begin();
      for (unsigned int i = 0; fit != bounds.end(); ++fit, ++i)
        {
        this->ROI[i] = (float)atof((*fit).c_str());
        }
      }
    else
      {
      PointListType seeds = this->GetSeeds();
      seeds[0];
      for (int i = 0; i < 3; i++)
        {
        this->ROI[2*i] = seeds[0].GetPosition()[i] - this->GetValueAsFloat("MaximumRadius");
        this->ROI[2*i+1] = seeds[0].GetPosition()[i] + this->GetValueAsFloat("MaximumRadius");
        }
      }
    return this->ROI;
    }
  
  itk::FixedArray< double, 3 > GetSigmas()
    {
    itk::FixedArray< double, 3 > s;
    std::list< std::string > bounds = this->GetValueAsList("ROI");
    std::list< std::string >::const_iterator fit = bounds.begin();
    for (unsigned int i = 0; fit != bounds.end(); ++fit, ++i)
      {
      s[i] = (double)atof((*fit).c_str());
      }
    return s;
    }

  PointListType GetSeeds()
    {
    typedef std::vector<std::vector<double> >  TempStorageType;
    TempStorageType tempPoints;
    std::vector<double> tempPoint;
    std::string seedFileName = this->GetValueAsString("SeedFile");
    std::ifstream inFile(seedFileName.c_str());
    std::string line;
    while( std::getline( inFile, line ) )
      {
      std::vector<std::string> strs;
      boost::split(strs, line, boost::is_any_of(","));
      tempPoint = std::vector<double>();
      for( unsigned int i = 0; i < strs.size(); ++i )
        {
        tempPoint.push_back(boost::lexical_cast<double>(strs[i]));
        }
      tempPoints.push_back(tempPoint);
      }
    PointListType seeds(tempPoints.size());
    for( unsigned int i = 0; i < tempPoints.size(); ++i )
      {
      seeds[i].SetPosition( tempPoints[i][0], tempPoints[i][1],
                            tempPoints[i][2] );
      }
    return seeds;
    }


protected:

  void AddArgument( std::string name, 
                    bool required,
                    std::string description,
                    TypeEnumType type = MetaCommand::STRING,
                    std::string defVal = "" )
    {
    this->SetOption(name,name,required,description);
    this->SetOptionLongTag(name,name);
    this->AddOptionField(name,name,type,true, defVal);
    } 

  double ROI[6];
};

#endif
