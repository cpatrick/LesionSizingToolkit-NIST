/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkImageToAIMXMLFilter.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageToAIMXMLFilter_h
#define __itkImageToAIMXMLFilter_h

#include "itkProcessObject.h"
#include "itkExtractImageFilter.h"
#include "itkContourExtractor2DImageFilter.h"

#include <map>
#include <vector>

namespace itk
{
  
/** \class ImageToAIMXMLFilter
 * \brief Converts an ITK image into an AIM XML file stored in a buffer.
 *
 *  This class takes takes the edges of a label image and exports it as the
 *  slice-level contours required by AIM.
 * 
 * \ingroup   ImageFilters     
 */
template <class TInputImage, class TReferenceImage >
class ITK_EXPORT ImageToAIMXMLFilter : public ProcessObject
{
public:
  /// Standard class typedefs.
  typedef ImageToAIMXMLFilter       Self;
  typedef ProcessObject             Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /// Method for creation through the object factory.
  itkNewMacro(Self);
  
  /// Run-time type information (and related methods).
  itkTypeMacro(ImageToAIMXMLFilter, ProcessObject);

  /// Constants
  static const int                                             Dimension = 2;

  /// Typedefs
  typedef TInputImage                                          InputImageType;
  typedef TReferenceImage                                      ReferenceImageType;
  typedef typename InputImageType::PixelType                   PixelType;
  typedef unsigned char                                        BinaryPixelType;
  typedef itk::Image<PixelType, Dimension>                     SliceType;
  typedef itk::Image<BinaryPixelType, Dimension>               ThresholdedSliceType;
  typedef typename InputImageType::ConstPointer                InputImagePointer;
  typedef typename ReferenceImageType::ConstPointer            ReferenceImagePointer;
  typedef std::vector<std::string>                             UIDContainerType;
  typedef itk::ExtractImageFilter<InputImageType, SliceType>   ExtractFilterType;
  typedef itk::ContourExtractor2DImageFilter<SliceType>        ContourFilterType;
  typedef std::vector<typename ReferenceImageType::IndexType>  ContourPointVectorType;
  typedef typename ReferenceImageType::IndexType::IndexValueType
                                                               ReferenceIndexValueType;
  typedef std::multimap<ReferenceIndexValueType,ContourPointVectorType>
                                                               ContourContainerType;

  /// Set the input in the form of an image (segmentation)
  void SetInput( const InputImageType * );
  
  /// Set the reference image for coordinate transformation
  void SetReference( const ReferenceImageType * );

  /// Filters the image into the AIM XML format
  void Update();

  /// Setup strings about the reference image from its DICOM headers
  itkSetStringMacro(CurrentUID);
  itkSetStringMacro(PatientName);
  itkSetStringMacro(PatientId);
  itkSetStringMacro(PatientSex);
  itkSetStringMacro(StudyInstanceUID);
  itkSetStringMacro(SeriesInstanceUID);

  /// Set the threshold at which to generate the contours on the input
  /// image
  itkSetMacro(ContourThreshold, PixelType);

  itkGetStringMacro(Output);

  /// Set the SOP Class UIDs that correspond to each slice of the reference
  /// image respectively
  void SetSOPClassUIDs( const UIDContainerType& uids );

  /// Set the SOP Instance UIDs that correspond to each slice of the reference
  /// image respectively
  void SetSOPInstanceUIDs( const UIDContainerType& uids );
  
protected:
  ImageToAIMXMLFilter(); 
  virtual ~ImageToAIMXMLFilter(); 

  /// Helper function for turning the vertex list into a contour to be inserted in
  /// the multimap m_Contours for extraction and output later
  void AddContourToVector( typename ContourFilterType::VertexListConstPointer verts,
                           unsigned int sliceNum,
                           ContourPointVectorType& vect,
                           ReferenceIndexValueType& index ) const;

  /// Helper function for iterating through the level-set segmentation passed
  /// to the filter and populating the multi-map of contour vectors
  void PopulateContourContainer();

  /// Helper function for generating the xml from the populated contour
  /// container and string metadata passed to the filter
  void GenerateXMLFromContours();

  /// Generate the current system time in the format required by AIM:
  /// YYYY-MM-DDTHH:MM:SS (the T is literal)
  std::string GetCurrentTimeInAIMFormat() const;

private:
  ImageToAIMXMLFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  InputImagePointer     m_InputImage;
  ReferenceImagePointer m_ReferenceImage;
  PixelType             m_ContourThreshold;
  std::string           m_Output;
  std::string           m_CurrentUID;
  std::string           m_PatientName;
  std::string           m_PatientId;
  std::string           m_PatientSex;
  std::string           m_StudyInstanceUID;
  std::string           m_SeriesInstanceUID;
  UIDContainerType      m_SOPClassUIDs;
  UIDContainerType      m_SOPInstanceUIDs;
  ContourContainerType  m_Contours;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToAIMXMLFilter.txx"
#endif

#endif
