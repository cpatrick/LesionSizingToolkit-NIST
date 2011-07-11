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
  /** Standard class typedefs. */
  typedef ImageToAIMXMLFilter       Self;
  typedef ProcessObject             Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToAIMXMLFilter, ProcessObject);

  /** Some typedefs. */
  static const int                                             Dimension = 2;
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

  /** Set the input in the form of an itk::Image */
  void SetInput( const InputImageType * );
  
  /** Set the reference image for coordinate transformation */
  void SetReference( const ReferenceImageType * );

  /** This filters the image into the AIM XML format */
  void Update();

  itkSetStringMacro(CurrentUID);
  itkSetStringMacro(PatientName);
  itkSetStringMacro(PatientId);
  itkSetStringMacro(PatientSex);
  itkSetStringMacro(StudyInstanceUID);
  itkSetStringMacro(SeriesInstanceUID);

  itkGetStringMacro(Output);

  void SetSOPClassUIDs( const UIDContainerType& uids );
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

private:
  ImageToAIMXMLFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  InputImagePointer     m_InputImage;
  ReferenceImagePointer m_ReferenceImage;
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
