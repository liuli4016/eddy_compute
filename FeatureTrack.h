#ifndef FEATURE_TRACK_H
#define FEATURE_TRACK_H

//#include "FeatureTrackUtil.h"

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <cassert>


using namespace std;


template <class T>
string TtoS(T arg)
{
  ostringstream os(ostringstream::out);
  os<<arg;
  return os.str(); // send to the ostringstream
}

template <class T>
T StoT(string arg)
{
  T a;

  istringstream(arg)>>a;
  return a;
}




struct Consts 
{
  static int const MAXSPLIT;
  static float const RADIUS;
  static float const FRAMEDIST;
  static float const NODEDIST;
  static long const MAXOBJS;
  static int const MAXLEN;
  static int const TRAKTABLE_MAXLEN; // important, for large dataset, this const might need to be changed. Otherwise, there will be segfault. Used in TrackObjects(),
  static float const DEFAULT_TOLERANCE;
  static bool const YES;
  static bool const NO;
  static int const DEFAULT_TIMESTEP_NUM;
  static bool const WriteTrackInfo;
/* If this value is set to 0 the TrackInfo file which contains
the tracking information table will not be written */
/* Set this value to 1 only if the entire set of timestep is being tracked */
/* Set the animate dual option to one step when feature tracking */
};



 

enum fieldtype { _UNIFORM_, _RECTILINEAR_,                   _UNSTRUCT_ };
enum DIRECTION {_FORWARD_, _BACKWARD_};
enum ROW_COL {_ROW_, _COL_};
enum SPLIT_MERGE {_SPLIT_, _MERGE_};
enum NEW_DISSIPATE {_NEW_, _DISSPATE_};
enum FILETYPE { _BASE_,_POLY_,_ATTR_,_UOCD_};
enum ATTRTYPE {_OBJNUM_,_INTEG_,_VOLUME_,_CENTX_,_CENTY_,_CENTZ_};

struct FT_NODE;
struct TrackObject;

struct TrackObject
{
     float time;
     int   volume;
     float mass;
     float xbar, ybar,zbar;
     float Ixx, Iyy, Izz, Ixy, Iyz, Izx; 
     float maximum;
     int  max_x, max_y, max_z;       
     bool existent;
     bool ismerge;
     bool issplit;
     bool isnew;
     bool flag;
     vector<int> parents; // seq_wr+no_init
     vector<int> children; // seq_wr+no_init 
     TrackObject (): volume(-1), mass(-1), existent(Consts::NO), isnew(Consts::NO), ismerge(Consts::NO), issplit(Consts::NO),flag(Consts::NO) 
     {
     }
     ~TrackObject()
     {
         vector<int>().swap(parents); 
         vector<int>().swap(children); 
     } 
};

struct ucdNode
{
     int ObjID;
     int NodeID;    
     inline friend int operator< (ucdNode t1, ucdNode t2) 
     {
         return (t1.NodeID<t2.NodeID?1:0);
     }
     inline friend ostream& operator<< (ostream& os,const ucdNode& t)
     {
          return os<<"("<<t.ObjID<<","<<t.NodeID<<")"<<" ";
     }
};


struct ObjectExtends
{
     long ObjVol; // the volume of the object
     double l_x; // lower left cordner coordinate of the surrounding boxof the object
     float l_y;
     float l_z;
     float u_x;// upper right corner coordinate
     float u_y;
     float u_z;
     float c_x; // centroid of the object
     float c_y;
     float c_z;
     float  minVal;
     float  maxVal;
    

};

struct Frame
{
     vector<ObjectExtends> objVols;  
  // [numObjs], noext+seq_wr+init, 
  //allocated in ReadTrak() with resize(), do not use push_back()!
  // objVols.size() is the number of objects in this frame
     vector<ucdNode> nodes; 

  // [numNodes], noext+rand_wr, 
  //allocated in ReadTrak() with resize(), do not use push_back()!
  // nodes.size() is the total volume of the objects in this frame.
  // Take easy, the unclassified nodes is good enough for overlap test.
     Frame()
     {
     }
     ~Frame()
     {
          vector<ObjectExtends>().swap(objVols);
          vector<ucdNode>().swap(nodes);
     }
};

struct FT_NODE
{ // means one object in one frame
     int Color[3];    // for polyfile coloring
     int FrameInd;  // starts from 0
     vector<int> Suc;//ext+rand_wr, extend with resize(n,-1)  
  //Suc should like "1 2 6 3 7 -1 -1 -1 -1 -1"
  //Cannot use size() to get the number of the valid elements, instead, use getnumSuc()
     vector<int> Pre; // ext+rand_wr, extend with resize(n,-1)  
     int numNodes;  // volume 
     int PredMaxnumNodes; // history max volume,  this is used to colorize
  //PreMaxnumNodes is necessary, it assures that the current object always be colored the same to the larger object. 
     int getnumSuc() const 
     {
          return std::count_if(Suc.begin(), Suc.end(), bind2nd(not_equal_to<int>(),-1));
     }
  
     int getnumPre() const 
     {
          return std::count_if(Pre.begin(), Pre.end(), bind2nd(not_equal_to<int>(),-1));
     }

     FT_NODE():FrameInd(-1),PredMaxnumNodes(0), Suc(vector<int>(Consts::MAXSPLIT,-1)), Pre(vector<int>(Consts::MAXSPLIT,-1))
     {
          Color[0]=Color[1]=Color[2]=-1;
     };
     ~FT_NODE() 
     {
          vector<int>().swap(Suc); // free all memory
          vector<int>().swap(Pre); // free all memory
     }
};

struct FRAME
{
     vector<FT_NODE> NodeArray; 
  //no_ext+seq_wr, actually seq_wr/rand_wr doesn't matter.  init in CopyTmpFrToList()
     int  index;
     inline FRAME():index(-1)
     {
     }
     inline FRAME(int const n):index(n)
     {
     }
     ~FRAME()
     {
          vector<FT_NODE>().swap(NodeArray);
     }
};

/////// For Colorize in viewDag.c/////////
struct TMPNODE
{
     int ObjInd;  // start from 0
     vector<int> Suc; 
  //ext+rand_wr, extend with resize(n,-1)  
     vector<int> Pre;
  //ext+rand_wr, extend with resize(n,-1) 

     TMPNODE():ObjInd(-1),Pre(vector<int>(Consts::MAXSPLIT,-1)),Suc(vector<int>(Consts::MAXSPLIT,-1)) 
     {
     }
     ~TMPNODE()
     {
          vector<int>().swap(Suc);
          vector<int>().swap(Pre);
     }
   
     int getnumSuc() const
     {
         return std::count_if(Suc.begin(),Suc.end(),bind2nd(not_equal_to<int>(),-1));
     }
 
     int getnumPre() const
     {
          return std::count_if(Pre.begin(),Pre.end(),bind2nd(not_equal_to<int>(),-1));
     }
};

struct TMPFR
{
     vector<TMPNODE> Obj; // ext+rand_wr 
  // NOTE: cannot use size() to get the number of valid elements
  // because there are some invalid ones.
  
     inline int getnumObjs()  const 
     {
          int i(0);
   
          for(vector<TMPNODE>::const_iterator it=Obj.begin();it!=Obj.end();++it)
          if(it->ObjInd!=-1) i++;
          return i;
     }
    
     TMPFR():Obj(vector<TMPNODE>(Consts::MAXOBJS))
     {
     }; 
  
     ~TMPFR()
     {
          vector<TMPNODE>().swap(Obj);
     }
};


template <class T>
void bound_check(vector<T> & v, int ind, int incre)
{
     int oldsize(v.size());
     int real_incre;
     if(v.size()>ind)
         return; // does't overflow
  
     if(ind==-1) //means directly increase the array by incre.
          real_incre=incre;
     else
          real_incre = incre > ind-v.size()+1 ? incre: ind-v.size()+1;
     int newsize(oldsize+real_incre);
     v.resize(newsize);
}


template <class T>
void bound_check(vector<T> & v, int ind, int incre, T initial_value)
{
     int oldsize(v.size());
     bound_check<T>(v,ind,incre);
     for(int i(oldsize);i<v.size();++i)  // now v.size is the new size
          v[i]=initial_value;
}

inline void objs_bcheck(int const step, int const obj,vector<vector<TrackObject> > & objs)
     // ad hoc for global variable objs
{
     bound_check<vector<TrackObject> >(objs,step,Consts::DEFAULT_TIMESTEP_NUM);
     bound_check<TrackObject> (objs[step], obj, Consts::MAXOBJS); // at
}

extern bool ReadList(char*  ListFileName,vector<string>& time_file, FILETYPE const which);
extern bool ReadTrak(string const baseName,Frame& frm);
extern bool ReadOct(string baseName,Frame& frm,int step,  vector<vector<TrackObject> > &objs);
extern int  createColormapFiles(string colormapfile,string label );
extern bool TrackObjects(string const basename,int const step,int curtime,Frame& f1,Frame& f2, vector<string>& time_polyfile,vector<vector<TrackObject> > & objs) ;
extern int  ReadPacketFile(string packetFName, vector<string> &previousFramePackInfo);




















extern void OverlapTest(Frame& t1,Frame& t2,vector<vector<int> > &OverlapTable,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2);
extern int ComputeScore(DIRECTION const direct,Frame& t1,Frame& t2,vector<vector<int> >& ScoreBoard,vector<vector<int> >& OverlapTable);
extern void FindOverlap(int const obj,vector<int>& Overlaps, DIRECTION const direct,Frame& t1,Frame& t2,vector<vector<int> > &OverlapTable );
extern void GenCombination(vector<int>& comb,int i);
extern int Intersect(vector<int> const & Comb,int const NumOvlp,vector<int> const & Overlaps,int const obj, DIRECTION const direct,vector<vector<int> > &OverlapTable);
extern float GeomMean(vector<int> const & Comb,int const NumOvlp,vector<int> const & Overlaps,int const obj, DIRECTION const direct, Frame& t1, Frame& t2);
extern bool RewritePolyFiles(string const trakTable, int const step,vector<string>& time_polyfile);
extern bool ReadDagFile(string const Filename,int const step, vector<FRAME>& FrameList);
extern void ParseTrackInfo(string const buffer,TMPFR &tmpfr0,TMPFR &tmpfr1);
extern bool AddNodeToTmpFr(int const which,int const obj,vector<int> & tmpObj,TMPFR &tmpfr0,TMPFR &tmpfr1);
extern void CopyTmpFrToList(int const which,TMPFR &tmpfr,vector<FRAME>& FrameList);
extern bool ReadNumNodes(vector<FRAME>& FrameList,int step,vector<string>& time_polyfile);
extern string path_core(string const s);
extern bool Read1stPolyFile(string const polyfile,vector<FRAME>& FrameList);
extern void Colorize(vector<FRAME>& FrameList);
extern void ColorizeSuc(vector<FT_NODE>::iterator node, vector<FRAME>::iterator nextfr);
extern bool UniformColorPoly(int const step,vector<FRAME>& FrameList,vector<string>& time_polyfile);
extern void TrackSplit_Merge(int const step, SPLIT_MERGE const sm,Frame& t1,Frame &t2, FILE* outfile,vector<vector<TrackObject> >&  objs,vector<int>& tag1,vector<int>& tag2,vector<vector<int> >& ScoreBoard, double &tolerance);
extern void TrackContinuous(int step,Frame& t1,Frame &t2,FILE* outfile,vector<vector<TrackObject> >&  objs,vector<int>& tag1,vector<int>& tag2,vector<vector<int> >& ScoreBoard, double& tolerance);
extern void TrackNew_Dissipate(int const step, NEW_DISSIPATE const nd,Frame& t1,Frame &t2,FILE* outfile,vector<vector<TrackObject> >&  objs,vector<int>& tag1,vector<int>& tag2);
extern void NumNonZeros(int const row,vector<int> & list, ROW_COL const rc,Frame& t1,Frame &t2 ,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2);
extern int GetSameScore(int const row,int const ind,int const nlist,vector<int> & list,vector<int> & split, ROW_COL const rc,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2, double &tolerance);
extern bool IsMax(int const row,int const col, ROW_COL const rc,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2,Frame& t1, Frame& t2, double &tolerance);
extern void GetInts(string const buffer, vector<int> & line);



#endif
