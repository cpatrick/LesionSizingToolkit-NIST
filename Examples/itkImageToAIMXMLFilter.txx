/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkImageToAIMXMLFilter.txx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkImageToAIMXMLFilter_txx
#define _itkImageToAIMXMLFilter_txx

#include "itkNumericTraits.h"

#include "itkImageToAIMXMLFilter.h"

namespace itk
{

/**
 * Constructor
 */
template <class TInputImage, class TReferenceImage>
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::ImageToAIMXMLFilter()
  :  m_InputImage(NULL),
     m_ReferenceImage(NULL),
     m_Output(""),
     m_CurrentUID(""),
     m_PatientName(""),
     m_PatientId(""),
     m_PatientSex(""),
     m_StudyInstanceUID(""),
     m_SeriesInstanceUID(""),
     m_SOPClassUIDs(),
     m_SOPInstanceUIDs(),
     m_Contours()
{
}

/**
 * Destructor
 */
template <class TInputImage, class TReferenceImage>
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::~ImageToAIMXMLFilter()
{
}

/**
 * Set an itk::Image as input 
 */
template <class TInputImage, class TReferenceImage>
void
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::SetInput( const InputImageType * inputImage )
{
  m_InputImage = inputImage;
}

/**
 * Set an itk::Image as input 
 */
template <class TInputImage, class TReferenceImage>
void
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::SetReference( const ReferenceImageType * referenceImage )
{
  m_ReferenceImage = referenceImage;
}

template <class TInputImage, class TReferenceImage>
void
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::SetSOPClassUIDs( const UIDContainerType& uids )
{
  m_SOPClassUIDs = uids;
}

template <class TInputImage, class TReferenceImage>
void
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::SetSOPInstanceUIDs( const UIDContainerType& uids )
{
  m_SOPInstanceUIDs = uids;
}
template <class TInputImage, class TReferenceImage>
void 
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::AddContourToVector( typename ContourFilterType::VertexListConstPointer verts,
                      unsigned int sliceNum,
                      ContourPointVectorType& vect,
                      ReferenceIndexValueType& index ) const
{
  typename ReferenceImageType::IndexType referenceIndex;
  for( unsigned int vertexNum = 0; vertexNum < verts->Size(); ++vertexNum )
    {
    std::cout << "Vertex: " << vertexNum << std::endl;
    typename ContourFilterType::VertexType vertex;
    typename InputImageType::IndexType inputIndex;
    typename InputImageType::PointType point;
    vertex = verts->ElementAt(vertexNum);
    inputIndex[0] = vertex[0];
    inputIndex[1] = vertex[1];
    inputIndex[2] = sliceNum;
    m_InputImage->TransformIndexToPhysicalPoint( inputIndex, point );
    m_ReferenceImage->TransformPhysicalPointToIndex( point, referenceIndex);
    for( unsigned int i = 0; i < 3; ++i )
      {
      std::cout << referenceIndex[i] << " ";
      }
    std::cout << std::endl;
    vect.push_back(referenceIndex);
    }
  index = referenceIndex[2];
}

/**
 * Delegate the Update to the importer
 */
template <class TInputImage, class TReferenceImage>
void
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::Update()
{
  typename ExtractFilterType::Pointer extractor = ExtractFilterType::New();
  extractor->SetInput( m_InputImage );
  typename InputImageType::RegionType inputRegion =
    m_InputImage->GetLargestPossibleRegion();
  typename InputImageType::SizeType inputSize = inputRegion.GetSize();
  typename InputImageType::IndexType inputIndex = inputRegion.GetIndex();
  typename InputImageType::SizeType desiredSize = inputSize;
  desiredSize[2] = 0;
  m_InputImage->Print( std::cout );
  std::cout << "**iterate" << std::endl;
  for( unsigned int sliceNum = inputIndex[2]; 
       sliceNum < inputSize[2]; 
       ++sliceNum )
    {
    std::cout << "Iteration: " << sliceNum << std::endl;
    typename InputImageType::RegionType desiredRegion;
    typename InputImageType::IndexType desiredIndex = inputIndex;
    typename ContourFilterType::Pointer contourer = ContourFilterType::New();
    inputIndex[2] = sliceNum;
    desiredRegion.SetSize(  desiredSize  );
    desiredRegion.SetIndex( desiredIndex );
    extractor->SetExtractionRegion( desiredRegion );
    extractor->SetDirectionCollapseToIdentity();
    contourer->SetInput(extractor->GetOutput());
    contourer->SetContourValue(-0.5);
    contourer->Update();
    for( unsigned int contourNum = 0;
         contourNum < contourer->GetNumberOfOutputs();
         ++contourNum )
      {
      std::cout << "Contour: " << contourNum << std::endl;
      typename ContourFilterType::VertexListConstPointer verts;
      verts = contourer->GetOutput(contourNum)->GetVertexList();
      ContourPointVectorType contourVector;
      IndexValueType referenceSlice;
      this->AddContourToVector( verts, sliceNum, contourVector, referenceSlice );
      typename ContourContainerType::value_type pairForInsertion(referenceSlice,
                                                                 contourVector);
      m_Contours.insert(pairForInsertion);
      }
    }
  
  m_Output = "<test></test>";
}

} // end namespace itk

#endif
