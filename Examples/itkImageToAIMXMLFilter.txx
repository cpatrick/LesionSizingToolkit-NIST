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
     m_SOPInstanceUIDs()
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
  for( int i = inputIndex[2]; i < inputSize[2]; ++i )
    {
    typename InputImageType::RegionType desiredRegion;
    typename InputImageType::IndexType desiredIndex = inputIndex;
    inputIndex[2] = i;
    desiredRegion.SetSize(  desiredSize  );
    desiredRegion.SetIndex( desiredIndex );
    extractor->SetExtractionRegion( desiredRegion );
    extractor->SetDirectionCollapseToIdentity();
    typename ThresholdFilterType::Pointer thresholder = ThresholdFilterType::New();
    thresholder->SetInput(extractor->GetOutput());
    thresholder->SetLowerThreshold(-0.5);
    thresholder->SetUpperThreshold(-0.5);
    thresholder->SetInsideValue(255);
    thresholder->SetOutsideValue(0);
    thresholder->Update();
    std::cout << i << std::endl;
    }
  
  m_Output = "<test></test>";
}

} // end namespace itk

#endif
