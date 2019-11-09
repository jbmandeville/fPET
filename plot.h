#ifndef plotData_H
#define plotData_H

#include <QMainWindow>
#include <QWidget>
#include <QVector>
#include "qcustomplot.h"
#include "io.h"

QT_BEGIN_NAMESPACE
class QAction;
class QGroupBox;
class QLabel;
class QMenu;
class QMenuBar;
class QPushButton;
class QFileDialog;
class QListWidget;
class QRadioButton;
class QString;
class QCPGraph;
QT_END_NAMESPACE

struct plotCurve
{
    dVector xData;        // original (non-concatenated x values); this should be removed from here and put into graphWindow
    dVector yData;        // data points
    bVector ignoreRescale;// if true, ignore the point when auto-scaling
    dVector yError;       // error bars (0=none, 1=bars, 2=envelope)
    dVector xPlot;        // potentially concatenated across files
    bool visible=true;    // allows curves to be turned off (e.g., single run mode)
    bool enabled=true;    // allows curves to be always invisible (but potentially available for export)
    bool histogram=false;
    bool exportCurve=true; // export curve to file? (e.g., ROIs)
    int iFile=0;        // keep track of which file (generally multiple files on graph for time series; concatenated or parallel)
    int iCurve;         // current curve number (index into _listOfCurves)
    double scaleFactor=1;
    bool isCurrentCurve=true; // for interactive cooperation between plotting and image display (curve = file, time point = tracer position)
    QString legend;
    QColor color=Qt::black;
    int lineThickness=1;
    int pointSize=5;
    QCPScatterStyle::ScatterShape pointStyle;
    int errorBars=0;      // 0=none, 1=bars, 2=envelop
};

class plotData : public QMainWindow
{
    Q_OBJECT
signals:
    void changedPointFromGraph(int plotID, int iFile, int iTime);

private:
    QCustomPlot *_qcplot;              // the QCustomPlot plotting class object
    QWidget *_parentPage;
    QStatusBar *_statusBar=nullptr;
    // Plotting mode and style
    int _plotID;                       // ID different plot surfacers (only the time plot gets some connections)
    bool _singleROIMode=true;          // if true, show only 1 ROI and potentially also the fit
    int _iFitType=0;                   // different fit type that indicate different plotting styles (data/fit, PET RTM with BP, .. etc)
    bool _autoScale=true;
    QCPRange _autoScaleXRange={0.,0.};
    bool _pressToMoveTracer=false;
    QCPItemTracer *_positionTracer;
    int _iFilePosition=0;
    int _iFilePositionLast=0;
    int _iTimePosition=0;
    int _iTimePositionLast=0;
    bool _concatenateRuns;

    bool _cursorOnPlot=false;

    int whichConcatenatedTracerFile(double x);
    int whichTracerTimePoint(int iFile);
    void interpretMousePosition(QMouseEvent *event);
    int getDataCurveIndex(int iFile);
    void concatenateRuns();

public:
    plotData(int plotID);
//    virtual ~plotData() {};
    QVector<plotCurve> _listOfCurves;
    // total number of curves = _listOfCurves.size() = _nFiles * _nCurvesPerFile + numberExtra
    int _nFiles=1;
    int _nCurvesPerFile=3;  // often data + fit + weights
    double _yAxis2Ratio;    // ratio of yAxis2 max to yAxis1 max

//    void updateListOfCurvesFromListOfFiles(QVector<ROIFile> ListOfFiles);
//    void updateListOfCurvesForTimeModel(ROIFile *CurrentROIFile, int iFile);
    void plotDataAndFit(bool newData);
    void setPlotStyle(int iStyle);    // inline getters
    void setPositionTracer(int iFile, double xPosition);
    void setPositionTracer(int iFile, int iTime);
    double setXPosition(int iFile, int iTime);
    void reScaleAxes(QCPRange *xRange, QCPRange *yRange);
    inline void setTracerDisplay(bool state) {_positionTracer->setVisible(state);}
    inline void setTracerXPosition(double xPos) {_positionTracer->setGraphKey(xPos);}
    void writeGraph(bool newGraph, QString fileName, QString regionName);

    inline QCustomPlot *getPlotSurface() const {return _qcplot;}
    inline int getFitType() const {return _iFitType;}
    inline bool getSingleROIMode() const {return _singleROIMode;}
    inline double getXPosition(int iCurve, int iTime) {return _listOfCurves[iCurve].xPlot[iTime];}

    // setting up a series of curves
    void init();
    void addCurve(int iFile, QString legend);
    void conclude(int iCurrentFile, bool singleRun);
    void setData(dVector xData, dVector yData);
    void setData(dVector xData, dVector yData, dVector yError);
    void setData(dVector xData, dVector yData, bVector ignoreRescale);
    void setData(dVector xData, dVector yData, dVector yError, bVector ignoreRescale);
    inline void setLineThickness(int thickness) {_listOfCurves.last().lineThickness=thickness;}
    inline void setPointSize(int pointSize) {_listOfCurves.last().pointSize=pointSize;
                                            if (_listOfCurves.last().pointStyle == QCPScatterStyle::ssNone)
                                                _listOfCurves.last().pointStyle = QCPScatterStyle::ssDisc;}
    inline void setPointStyle(QCPScatterStyle::ScatterShape pointStyle) {_listOfCurves.last().pointStyle=pointStyle;}
    inline void setErrorBars(int errorBars) {_listOfCurves.last().errorBars=errorBars;}
    inline void setColor(QColor color) {_listOfCurves.last().color=color;}
    inline void setVisible(bool visible) {_listOfCurves.last().visible=visible;}
    inline void setEnabled(bool enabled) {_listOfCurves.last().enabled=enabled;}
    inline void setHistogram(bool histogram) {_listOfCurves.last().histogram=histogram;}
    inline void setExport(bool exportCurve) {_listOfCurves.last().exportCurve=exportCurve;}
    inline void setCurrentCurve(bool isCurrentCurve) {_listOfCurves.last().isCurrentCurve=isCurrentCurve;}
    inline void setLegend(QString legend) {_listOfCurves.last().legend=legend;}
    inline void setScaleFactor(double scaleFactor) {_listOfCurves.last().scaleFactor=scaleFactor;}
    inline void setYAxisRatio(double ratio) {_yAxis2Ratio = ratio;}

    inline void setQCStatusBar(QStatusBar *statusBar) {_statusBar = statusBar;}
    inline void setAutoScale( bool value ) {_autoScale = value;}
    inline void setSingleROIMode( bool mode ) {_singleROIMode = mode;}
    inline void setVisibleXAxis(bool visible) {_qcplot->xAxis->setVisible(visible);}
    inline void setLabelXAxis(QString label) {_qcplot->xAxis->setLabel(label);}
    inline void setLabelYAxis(QString label) {_qcplot->yAxis->setLabel(label);}
    inline void setFitType(int iFit) {_iFitType = iFit;}
    inline void setMainPage(QWidget *main) {_parentPage = main;}
    inline void setAutoScaleXRange(QCPRange range) {_autoScaleXRange = range;}
    inline void setLegendOn(bool state) {_qcplot->legend->setVisible(state);}
    inline void setConcatenatedRuns(bool state) {_concatenateRuns = state; concatenateRuns();}

    // getters
    inline double getTracerFile() {return _iFilePosition;}
    inline double getTracerTimeIndex() {return _iTimePosition;}
    inline double getTracerXPosition() {return _positionTracer->graphKey();}
    inline QString getLabelXAxis() {return _qcplot->xAxis->label();}
    double getMaxYAbs(plotCurve *ptrCurve);
    inline void setXRange(QCPRange range) {_qcplot->xAxis->setRange(range); _qcplot->replot();}
    inline void setYRange(QCPRange range) {_qcplot->yAxis->setRange(range); _qcplot->replot();}
    inline bool isAutoScale() {return _autoScale;}
    inline plotCurve *getThisCurvePointer() {return &_listOfCurves.last();}
    inline int getNumberCurves() {return _listOfCurves.size();}
    inline dVector getXData(int iCurve) {return _listOfCurves[iCurve].xData;}
    inline dVector getYData(int iCurve) {return _listOfCurves[iCurve].yData;}

private slots:
    void axisLabelDoubleClick(QCPAxis *axis, QCPAxis::SelectablePart part);
    void legendDoubleClick(QCPLegend *legend, QCPAbstractLegendItem *item);
    void contextMenuRequest(QPoint pos);
    void makeSelectedGraphDotted();
    void makeSelectedGraphBiggerPoints();
    void makeSelectedGraphSmallerPoints();
    void removeSelectedGraph();
    void writeSelectedGraph();
    void moveLegend();
    void plotMousePress(QMouseEvent *event);
    void plotMouseMove(QMouseEvent *event);
    void writeOneGraph();
    inline void showLegend() {_qcplot->legend->setVisible(true); _qcplot->replot();}
    inline void hideLegend() {_qcplot->legend->setVisible(false); _qcplot->replot();}

    void lastTimePoint();
    void setAxis2Range(QCPRange Range);

    void keyboardPlus();
    void keyboardMinus();

public slots:
    void setXZoom();
    void setYZoom();
    void autoScale(bool state);
    void setSelectPoints();
    void changePoint(int plotID, int iFile, int iTime);

};

#endif // plotData_H
