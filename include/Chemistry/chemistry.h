/****************************************************************************
 *                                                                          *
 * Copyright (C) 2015 Neutrino International Inc.                           *
 *                                                                          *
 * Author : Brian Lin <lin.foxman@gmail.com>, Skype: wolfram_lin            *
 *                                                                          *
 ****************************************************************************/

#ifndef QT_CHEMISTRY_H
#define QT_CHEMISTRY_H

#include <QtManagers>

QT_BEGIN_NAMESPACE

#ifndef QT_STATIC
#    if defined(QT_BUILD_CHEMISTRY_LIB)
#      define Q_CHEMISTRY_EXPORT Q_DECL_EXPORT
#    else
#      define Q_CHEMISTRY_EXPORT Q_DECL_IMPORT
#    endif
#else
#    define Q_CHEMISTRY_EXPORT
#endif

namespace N
{

namespace ProteinBank
{

class Q_CHEMISTRY_EXPORT AtomicPoint ;
class Q_CHEMISTRY_EXPORT Atom        ;
class Q_CHEMISTRY_EXPORT Molecule    ;
#ifdef Q_OS_WIN
class Q_CHEMISTRY_EXPORT PdbViewer   ;
#endif

class Q_CHEMISTRY_EXPORT AtomicPoint
{
  public:

                  AtomicPoint (void) ;
                  AtomicPoint (float x,float y,float z) ;
    virtual      ~AtomicPoint (void) ;

    float         x           (void) const ;
    float         y           (void) const ;
    float         z           (void) const ;
    const float * values      (void) const ;

    void          setX        (float x) ;
    void          setY        (float y) ;
    void          setZ        (float z) ;

  protected:

    float Values [ 3 ] ;

  private:

};

class Q_CHEMISTRY_EXPORT Atom : public AtomicPoint
{
  public:

             Atom    (void) ;
    virtual ~Atom    (void) ;

    int      code    (void) const ;
    void     setCode (int code) ;
    void     setCode (const QString & strCode) ;

  protected:

    int Code ;

};

class Q_CHEMISTRY_EXPORT Molecule : QObject
{
  public:

    typedef QPair<int,int> Link ;

                 Molecule (QObject * parent = NULL) ;
    virtual     ~Molecule (void) ;

    void         load     (const QString    & filename) ;
    void         load     (const QByteArray & body    ) ;
    void         parse    (QTextStream      & stream  ) ;
    int          numAtoms (void) const ;
    int          numLinks (void) const ;
    const Atom & getAtom  (int id) const ;
    const Link & getLink  (int id) const ;

  protected:

    QVector<Atom> m_Atoms ;
    QVector<Link> m_Links ;

    void         addLink  (int id1,int id2) ;

  private:

  public slots:

  protected slots:

  private slots:

  signals:

};

#ifdef Q_OS_WIN

class Q_CHEMISTRY_EXPORT PdbViewer : public GLWidget
{
  Q_OBJECT
  public:

    static const float cpkTable[16][4] ;

    enum Style                {
      STYLE_CPK_SPHERES       ,
      STYLE_SPHERES_CYLINDERS ,
      STYLE_CYLINDERS         ,
      STYLE_WIRES           } ;

    struct Link         {
      int     atomId    ;
      float   height    ;
      GLfloat rotVec[3] ;
      float   angle     ;
    }                   ;

                 PdbViewer         (StandardConstructor) ;
    virtual     ~PdbViewer         (void) ;

    virtual QSize sizeHint         (void) const ;

    void         setMolecule       (const Molecule *mol);

  protected:

    GLUquadricObj  * quadric    ;
    QPoint           lastPos    ;
    const Molecule * molecule   ;

    virtual void initializeGL      (void) ;
    virtual void paintGL           (void) ;
    virtual void resizeGL          (int width,int height) ;

    virtual void focusInEvent      (QFocusEvent       * event) ;
    virtual void focusOutEvent     (QFocusEvent       * event) ;
    virtual void contextMenuEvent  (QContextMenuEvent * event) ;
    virtual void mousePressEvent   (QMouseEvent       * event) ;
    virtual void mouseMoveEvent    (QMouseEvent       * event) ;
    virtual void mouseReleaseEvent (QMouseEvent       * event) ;
    virtual void wheelEvent        (QWheelEvent       * event) ;
    virtual void keyPressEvent     (QKeyEvent         * event) ;
    virtual void resizeEvent       (QResizeEvent      * event) ;
    virtual void showEvent         (QShowEvent        * event) ;

  private:

    int           xRot       ;
    int           yRot       ;
    int           zRot       ;
    float         distance   ;
    float         xMin       ;
    float         xMax       ;
    float         yMin       ;
    float         yMax       ;
    float         zMin       ;
    float         zMax       ;
    bool          flagBox    ;
    bool          flagAxis   ;
    Style         flagStyle  ;
    QVector<Link> links      ;
    QColor        background ;

    float getRotationAngle         (const AtomicPoint & coord1  ,
                                    const AtomicPoint & coord2  ,
                                    float             * vrot  ) ;
    void  renderAtoms              (void) ;
    void  renderLinks              (void) ;
    void  renderWires              (void) ;
    void  renderBox                (void) ;
    void  renderAxis               (void) ;

  public slots:

    void showBox                   (bool enabled) ;
    void showAxis                  (bool enabled) ;
    void setStyle                  (Style estilo) ;

    virtual void Export            (void) ;
    virtual void startup           (void) ;
    virtual void Relocation        (void) ;

  protected slots:

    virtual bool Menu              (QPoint pos) ;
    virtual void distanceChanged   (double v) ;

  private slots:

  signals:

};

#endif

}

}

QT_END_NAMESPACE

#endif
