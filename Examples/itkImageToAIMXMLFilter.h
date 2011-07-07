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
#include "vtkImageImport.h"
#include "vtkImageData.h"

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
template <class TInputImage >
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
  typedef TInputImage                                 InputImageType;
  typedef typename    InputImageType::ConstPointer    InputImagePointer;
 

  const char* GetOutput() const;

  /** Set the input in the form of an itk::Image */
  void SetInput( const InputImageType * );

  /** This filters the image into the AIM XML format */
  void Update();
  
protected:
  ImageToAIMXMLFilter(); 
  virtual ~ImageToAIMXMLFilter(); 

private:
  ImageToAIMXMLFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  const char*       m_XMLBuffer;
  InputImagePointer m_InputImage;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToAIMXMLFilter.txx"
#endif

#endif
