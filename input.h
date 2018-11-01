#ifndef INPUT_H
#define INPUT_H
// ************************************************************************* //
//  File: input.h
// ************************************************************************* //
#include <string>
#include <vtkRectilinearGrid.h>
#include <vtkDataSet.h>
#include "nodeinfo.h"



#define METHOD_SUCCESS 1
#define METHOD_FAILURE 0
#define NOT_FOUND  0
#define YES_FOUND 1

const int _PRISM_=0;
const int _HEX_=1;
const int _QUAD_=2;

const int _NONUNIFORM_=1;
const int THRESH_PERCENT = 0;
const int INVALID_OBJ = -1;
const int INTERVAL = 1000;
const int ADAPTIVE = 1;
const int FIXED = 0;

#define  F_UNUSED  0x00000000
#define  F_THR     0x00000001
#define  F_USED    0x00000002
#define  F_UCD     0x00000004





struct stRGB
{
  int r, g, b;
};


struct stObject
{
  stCellIndex *cellPtr; 
  stObject    *next;
  stCellIndex *endPtr;
  long        objNum;
  long        numCell; 
  long        maxNode;
  long        minNode;
  float       objMax;
  float       objMin;
  float       minX;
  float       minY;
  float       minZ;
  float       maxX;
  float       maxY;
  float       maxZ;
  stRGB       color;
  long        volume;
  float       mass;
  float       massSqr;
  float       centroidX;
  float       centroidY;
  float       centroidZ;
  float       Ixx;
  float       Iyy;
  float       Izz;
  float       Ixy;
  float       Iyz;
  float       Izx;

  float       LowerLeft_Corner_x;
  float       LowerLeft_Corner_y;
  float       LowerLeft_Corner_z;
  float       UpperRight_Corner_x;
  float       UpperRight_Corner_y;
  float       UpperRight_Corner_z;
  long        Packet_ID;
  bool        packetFlag;
};







template <class CELLtype>
class InpObject 
{ 
 public:
  static int celltype;	
  char fileName[256];
  unsigned long numObj;
  unsigned long numCell;
  unsigned long numNodes;
  int  colorMethod;
  float minData;
  float maxData;
  float minData1;
  float maxData1;
  float minData2;
  float maxData2;  
  float thrVal;
  float thrVal1;
  float thrVal2;  
  float thresh_deltax;
  float thresh_deltay;
  float thresh_deltaz;
  int   thrType;
  float threshPercent_FromTop;
  float threshPercent_FromBottom;
  int nncomp;
  long MaxPacketnumber; //TotalPacket number in the current frame
  CELLtype  *pcell; // simple free
  stNodePos    *pnode; // complex free
  stNodeData   *pnodeData; //simple free
  stObject *pobject; // complex free

  stPacket *packets;    // contains the packets 

  InpObject() {pcell=0; pnode=0; pnodeData=0; pobject=0; packets=0;}
  ~InpObject(); 
  int freeCells(stCellIndex*);
  int freePacketLists(ObjIndex*);
};

//*** +++


template<class T>
InpObject<T> :: ~InpObject() 
{
     if(pcell)
          delete [] pcell;
     if(pnodeData)
          delete [] pnodeData;  
     if(packets)
          delete [] packets;      

     stObject *objPtr, *nextPtr; // freeObjects
     stCellIndex   *cellPtr;
     objPtr=pobject;
       

    while(objPtr)
    {
         nextPtr=objPtr->next;
  	 cellPtr=objPtr->cellPtr;
  	 if(cellPtr)
              freeCells(cellPtr);
  	 delete objPtr;
	 objPtr=nextPtr;  	
    }
    stCellIndex *list, *nextlist; // freePos
    register int i;
    for(i=0;i<numNodes;i++)
    {
         delete pnode[i].adjPosList;
         list=pnode[i].list;
         while(list)
         {
              nextlist=list->next;
              delete list;
              list = nextlist;
         }
    }
    delete [] pnode;

}



template<class T>
int InpObject<T>::freePacketLists(ObjIndex *ObjPtr )
{
     ObjIndex    *nextPtr;    
     while(ObjPtr)
     {
          nextPtr=ObjPtr->next;
          delete ObjPtr;
          ObjPtr=nextPtr;	
     } 
  return METHOD_SUCCESS; 
}



template<class T>
int InpObject<T>::freeCells( stCellIndex *cellPtr )
{
     stCellIndex    *nextPtr;    
     while(cellPtr)
     {
          nextPtr=cellPtr->next;
          delete cellPtr;
          cellPtr=nextPtr;	
     } 
  return METHOD_SUCCESS; 
}





//*** ---
template <class UCDtype>
class InpObject;



template <class T>
int AddObjToObjList( InpObject<T> *inPtr, stObject *obj, long vol,int minObjSize)
{
     stObject *ptr;
     
     if(vol<minObjSize)          /* filter test */
          return METHOD_SUCCESS; 
     else                        /* changed 11/97  */
     {   
          ptr = inPtr->pobject;
      
         if (ptr == NULL)          /* no prior object */
         {
             inPtr->pobject = obj;   /* link first object */
	     inPtr->numObj=0;
	     obj->objNum=0;
	     obj->next=NULL;
	     inPtr->numObj++;
       	     return METHOD_SUCCESS;
         }
      
         else                      /* append obj to linked list */
	 {
	      while( ptr->next != NULL )
	      {
	           ptr=ptr->next;
	      }
	      ptr->next=obj;
	      obj->next = NULL;
	      obj->objNum=inPtr->numObj;
	      inPtr->numObj++;  /* increment num of objects */
             return  METHOD_SUCCESS;
	 }
    }  
    return  METHOD_SUCCESS;
}



// ************************************************************************* //
//  END: inputdata.h
// ************************************************************************* //

#endif
