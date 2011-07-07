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
template <class TInputImage>
ImageToAIMXMLFilter<TInputImage>
::ImageToAIMXMLFilter()
{
}

/**
 * Destructor
 */
template <class TInputImage>
ImageToAIMXMLFilter<TInputImage>
::~ImageToAIMXMLFilter()
{
}

/**
 * Set an itk::Image as input 
 */
template <class TInputImage>
void
ImageToAIMXMLFilter<TInputImage>
::SetInput( const InputImageType * inputImage )
{
  m_InputImage = inputImage;
}

/**
 * Get the xml buffer of the output
 */
template <class TInputImage>
const char *
ImageToAIMXMLFilter<TInputImage>
::GetOutput() const
{
  return m_XMLBuffer;
}

/**
 * Delegate the Update to the importer
 */
template <class TInputImage>
void
ImageToAIMXMLFilter<TInputImage>
::Update()
{
  m_XMLBuffer = "<test></test>";
}

} // end namespace itk

#endif
