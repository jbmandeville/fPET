#ifndef IO_H
#define IO_H

#include <QString>
#include <QVector>
#include <QProcess>
#include <QtMath>

#define PI 3.141592654
#define PI_180 (PI/180.)

#define NCOLORS 14 // do not include undefined value
enum overlayColor
{   //                 R   G   B
    Color_Cyan,    //  0  255 255
    Color_Teal,    //  0  128 128
    Color_Purple,  // 128 128  0
    Color_Magenta, // 255  0  255
    Color_Pink,    // 255 192 203
    Color_Red,     // 255  0   0
    Color_Lime,    //  0  255  0
    Color_Green,   //  0  128  0
    Color_Sky,     //  0  191 255
    Color_Blue,    //  0   0  255
    Color_Yellow,  // 255 255  0
    Color_Orange,  // 255 165  0
    Color_Tan,     // 210 180 140
    Color_Brown,   // 139 69  19
    Color_Undefined // this must be last, so that the enum values match the combobox
};

enum fileStorageTypes
{
    BFLOAT    = 0,
    BSHORT    = 1,
    BLONG     = 2,
    UINT8     = 3,
    USHORT    = 4,
    BDOUBLE   = 5,
    NIFTI_NII = 6,
    NIFTI_IMG = 7
};

enum fileCategory
{
    category_none,
    category_timeSeries,
    category_anatomy,
    category_atlas,
    category_preProcessingSource,
    category_color,
    category_pValue,
    category_TValue,
    category_signalValue,
    category_BPValue,
    category_ADC,
    category_anisotropy,
    category_vector,
    fileCategory_MAX
};

struct IPoint2D {union {int x; int lower;}; union {int y; int upper;}; };
struct IPoint3D {int x, y, z;};
struct IPoint4D {int x, y, z, t;};
struct Point2D  {float x, y;};
struct Point3D  {float x, y, z;};
struct Point4D  {float x, y, z, t;};
struct dPoint2D {double x, y;};
struct dPoint3D {double x, y, z;};
struct dPoint4D {double x, y, z, t;};
struct Mat33 {double m[3][3];};
struct Mat44 {double m[4][4];};
struct fComplex {float real, imag;};
struct dComplex {double real, imag;};
typedef QVector<int>                iVector;
typedef QVector<bool>               bVector;
typedef QVector<QChar>              cVector;
typedef QVector<QString>            sVector;
typedef QVector<double>             dVector;
typedef QVector<float >             fVector;
typedef QVector<IPoint2D>          i2Vector;
typedef QVector<IPoint3D>          i3Vector;
typedef QVector<dPoint3D>          d3Vector;
typedef QVector<fComplex>          fcVector;
typedef QVector<dComplex>          dcVector;
typedef QVector<QVector<bool>>      bMatrix;
typedef QVector<QVector<QChar>>     cMatrix;
typedef QVector<QVector<int>>       iMatrix;
typedef QVector<QVector<float>>     fMatrix;
typedef QVector<QVector<double>>    dMatrix;
typedef QVector<QVector<fComplex>> fcMatrix;
typedef QVector<QVector<Point3D>>  f3Matrix;
typedef QVector<QVector<dPoint3D>> d3Matrix;
typedef QVector<QVector<IPoint2D>> i2Matrix;
typedef QVector<QVector<QVector<int>>>    iMatrix3;
typedef QVector<QVector<QVector<double>>> dMatrix3;
struct fVoxel {int x, y, z; float value; overlayColor color;};  // can hold a voxel location, value, and optionally a color for display
typedef QVector<fVoxel> VoxelSet;

struct OVLFile
{  // typical overlay-list format: [name] [path] [optional color]
    IPoint3D _dimensions={0,0,0};  // keeps image dimensions (not always available)
    QString name;                  // e.g. putamen
    QString path;                  // e.g. $TemplateDir/putamen.ovl
    overlayColor color = Color_Cyan;  // Attach a color to an overlay file
    int iAtlasRegionID;  // indicates region from atlas, with ID as given (e.g region 21); set to 0 for non-atlas regions
    VoxelSet voxelList;// vector of (x,y,z,value)
    bool loaded=false;
};

struct imageHeader
{ // all info to interpret a file should be stored in the header in order to facilitate voxel-reading in a different thread
    QString fileName;       // original file name (if read from file), including directory path
    IPoint4D dim={0,0,0,0}; // # voxels in x, y, z, t
    dPoint4D resolution={0.,0.,0.,0.};    // resolutions in x, y, z
    IPoint4D points={0,0,0,0};        // (nx, nx*ny, nx*ny*nz, nx*ny*nz*nt)
    dMatrix ijk_to_xyz;     // 4x4 transformation from data indices (i,j,k) to space (x,y,z)
    dMatrix xyz_to_ijk;     // 4x4 inverse transformation
    int byte_order;         // byte order 0 = Mac, 1 = PC
    bool offDiagonal=false ;// if true, the transformation matrix does not have parallel voxels & coordinates
    int fileType=NIFTI_NII; // e.g., bshort, bfloat, blong, nifti_nii, nifti_img
    int dataType=0;         // 0 = float, 1 = real-imaginary, 2 = magnitude-phase, 3 = vector
    int storageType=BFLOAT; // storage type: bshort, bfloat, blong
    double DOF=0.;          // degrees of freedom for statistics; default value 0
    int voxelOffset;        // offset to the first voxel data in the file
    float scaleSlope;       // scaling factor for data (data --> slope*data + offset); ignore if 0; default 0
    float scaleOffset;      // offset for data; default 0
};

struct voxelStack
{
    fMatrix  f1;      // 4-dimensional data vector; allocated as an array of pointers to volumes; dimensions [time][i3d]
    fcMatrix fc;      // 4-dimensional complex vector
    f3Matrix f3;      // 4-dimensional 3-vector
    iVector color;    // 3-dimensional color (for bitmaps/overlays)
    bool empty=false; // if true, all data is zero
};

// Define a structure for a stack of images
struct imageData
{
    voxelStack voxels;
    imageHeader hdr;

//    IPoint4D dim={0,0,0,0}; // # voxels in x, y, z, w=t
//    Point4D resolution;    // resolutions in x, y, z, t
//    Point4D origin;        // origins in x, y, z, t
//    Point3D direction;     // directions in x, y, z
//    IPoint4D points;       // (nx, nx*ny, nx*ny*nz, nx*ny*nz*nt)
//    Mat44 ijk_to_xyz;      // transformation from data indices (i,j,k) to space (x,y,z)
//    Mat44 xyz_to_ijk;      // inverse transformation
};

struct ROI_data
{
    QString fileName;// Vector of file names for each run in the ROI
    QString name;    // some identifier
    IPoint3D voxel;  // voxel for 1-point ROI
    double nspace;   // # points in space (non-integer due to potential non-binary voxel weights)
    double dt;       // time step
    double dof;
    dVector xTime;   // 1-dimensional data vector [itime]
    dVector ySignal; // 1-dimensional data vector [itime]
};

namespace utilIO
{
    template <class T> T* allocate_vector(int n)
    {
        T* a = new T[n];
        return a;
    }
    template <class T> void delete_vector(T *vec)
    {
        delete vec;
    }

    int machine_byte_order();
    void swap(short *ptr);
    void swap(unsigned short *ptr);
    void swap(float *ptr);
    void swap(int  *ptr);
    void swap_16( void *ptr );
    void swap_32( void *ptr );

    inline int index3d(IPoint3D dim, IPoint3D voxel) {return voxel.z*dim.x*dim.y + voxel.y*dim.x + voxel.x;}
    inline int index3d(IPoint3D dim, int iX, int iY, int iZ) {return iZ*dim.x*dim.y + iY*dim.x + iX;}

    int readOverlayFile(VoxelSet &voxelList, QString fileName, overlayColor color, IPoint3D dimSpace);

    void delayMS( int millisecondsToWait );

    int readTableFile(int iRun, QString fileName, QStringList &columnNames, dMatrix3 &table, QStringList &warningErrors);

}

namespace utilMath
{
#define SQR(a) ((a) == 0.0 ? 0.0 : (a)*(a))
#define NMAT 3

    void FFTW1D_Volume( IPoint3D dim, dcVector &volume, int iDimFFT, bool forward );

    Mat44 invertMat44(Mat44 matrix44 );
    bool dInvertSquareMatrix(dMatrix &dSquareMatrix );

    double deltaKernel(double x);
    double LanczosKernel(double x, double width);
    double GaussKernel(double x, double y, double fwhm);
    double inverse_Gauss(double x, double y, double fwhm);
    void swapX( dComplex *dcData, int lDim, int lAdd );
    bool ParabolicInterpolation(double *xParabola, double *yParabola, double &xMax, double &yMax);
    void fuzzyBinning(double value, double min, double max, int nBins, IPoint2D &iBin, dPoint2D &weightBin );

    void matrixProduct( dMatrix X1, dMatrix X2, dMatrix &X1X2 );
    void matrixProductTranspose1( dMatrix X1, dMatrix X2, dMatrix &X1X2 );
    void matrixProductTranspose2( dMatrix X1, dMatrix X2, dMatrix &X1X2 );
    void matrixPseudoInverse( dMatrix X_12, dMatrix &piX_21 );
    double traceOfSquareMatrix(dMatrix mat);

    void eigen_decomposition(double A[NMAT][NMAT], double V[NMAT][NMAT], double d[NMAT]);
    inline double hypot2(double x, double y) {return qSqrt(x*x+y*y);}
    void tred2(double V[NMAT][NMAT], double d[NMAT], double e[NMAT]);
    void tql2(double V[NMAT][NMAT], double d[NMAT], double e[NMAT]);

    double ibeta(double aa, double bb, double xx);
    double incbcf(double a, double b, double x);
    double incbd(double a, double b, double x);
    double pseries(double a, double b, double x);

    void topDownMergeSort(dVector &array);
    void topDownSplitMerge(int iBegin, int iEnd, dVector &array, dVector &arrayWorking);
    void topDownMerge(int iBegin, int iMiddle, int iEnd, dVector &array, dVector &arrayWorking);
    void copyArray(int iBegin, int iEnd, dVector &array, dVector &arrayWorking);

    void topDownMergeSort(dVector &array, iVector &indexArray);
    void topDownSplitMerge(int iBegin, int iEnd, dVector &array, dVector &arrayWorking, iVector &indexArray, iVector &indexWorking);
    void topDownMerge(int iBegin, int iMiddle, int iEnd, dVector &array, dVector &arrayWorking, iVector &indexArray, iVector &indexWorking);
    void copyArray(int iBegin, int iEnd, dVector &array, dVector &arrayWorking, iVector &indexArray, iVector &indexWorking);

    VoxelSet createVoxelList(int iThread, int nThreads, IPoint3D dim);
}
namespace utilString
{
    void errorMessage(QString errorText);
    int decodeSelectionList(QString list, iVector &includeVolume);
    int decodeSelectionList(QString list, bVector &includeVolume);
    QString recodeSelectionList(iVector includeVolume);
    QString recodeSelectionList(bVector includeVolume);
    void replaceEnvironmentVariables(QStringList &inputStringList);
    QString replaceEnvironmentVariables(QString inputString);
    QString insertEnvVariableTemplateDir(QString inputString);
    QString getFileNameExtension(QString fileName);
    QString getFileNameWithoutExtension(QString fileName);
    QString getFileNameWithoutDirectory(QString fileName);
    QString getFileNameBase(QString fileName);
    QString getDirectoryName(QString fileName);
    bool fileHasExtension(QString fileName);
}
#endif // IO_H
