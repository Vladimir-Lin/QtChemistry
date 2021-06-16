#include <chemistry.h>

#ifdef Q_OS_WIN

#define CX 0
#define CY 1
#define CZ 2

typedef struct                  {
  char         symbol [ 3 ]     ;
  int          covalentRadio    ;
  int          vanDerWaalsRadio ;
  int          cpkColor         ;
  const char * name             ;
} AtomInfo                      ;

extern AtomInfo Atoms [ ]       ;

float Distance(const float p1[],const float p2[])
{
  float dx = p1[CX] - p2[CX]                    ;
  float dy = p1[CY] - p2[CY]                    ;
  float dz = p1[CZ] - p2[CZ]                    ;
  return ::sqrt ( dx * dx + dy * dy + dz * dz ) ;
}

float vectorialModulus(const float * vec)
{
  return ::sqrt ( vec[CX] * vec[CX] + vec[CY] * vec[1] + vec[CZ] * vec[2] ) ;
}

float vectorialAngle(const float * A,const float * B)
{
  float prod = A[CX] * B[CX] + A[CY] * B[CY] + A[CZ] * B[CZ] ;
  float mod1 = vectorialModulus ( A )                        ;
  float mod2 = vectorialModulus ( B )                        ;
  const float PI = ::acos ( -1.0f )                          ;
  float rad = acos ( prod / ( mod1 * mod2 ) )                ;
  float deg = rad * 180 / PI                                 ;
  return deg                                                 ;
}

void vectorialProduct(const float * A,const float * B,float * prod)
{
  prod[CX] =     A[CY] * B[CZ] - A[CZ] * B[CY]   ;
  prod[CY] = - ( A[CX] * B[CZ] - A[CZ] * B[CX] ) ;
  prod[CZ] =     A[CX] * B[CY] - A[CY] * B[CX]   ;
}

#define ATOM_RADIO    0.5f
#define LINK_RADIO    0.18f

const float N::ProteinBank::PdbViewer::cpkTable[16][4] = {
    { 200/255.0f, 200/255.0f, 200/255.0f, 1.0f }, /* Light Grey */
    { 143/255.0f, 143/255.0f, 255/255.0f, 1.0f }, /* Sky Blue */
    { 240/255.0f,   0/255.0f,   0/255.0f, 1.0f }, /* Red */
    { 255/255.0f, 200/255.0f,   0/255.0f, 1.0f }, /* Yellow */
    { 255/255.0f, 255/255.0f, 255/255.0f, 1.0f }, /* White */
    { 255/255.0f, 192/255.0f, 203/255.0f, 1.0f }, /* Pink */
    { 218/255.0f, 165/255.0f,  32/255.0f, 1.0f }, /* Golden Rod */
    {   0/255.0f,   0/255.0f, 255/255.0f, 1.0f }, /* Blue */
    { 255/255.0f, 165/255.0f,   0/255.0f, 1.0f }, /* Orange */
    { 128/255.0f, 128/255.0f, 144/255.0f, 1.0f }, /* Dark Grey */
    { 165/255.0f,  42/255.0f,  42/255.0f, 1.0f }, /* Brown */
    { 160/255.0f,  32/255.0f, 240/255.0f, 1.0f }, /* Purple */
    { 255/255.0f,  20/255.0f, 147/255.0f, 1.0f }, /* Deep Pink */
    {   0/255.0f, 255/255.0f,   0/255.0f, 1.0f }, /* Green */
    { 178/255.0f,  34/255.0f,  34/255.0f, 1.0f }, /* Fire Brick */
    {  34/255.0f, 139/255.0f,  34/255.0f, 1.0f }, /* Forest Green */
};

N::ProteinBank::PdbViewer:: PdbViewer ( QWidget * parent , Plan * p      )
                          : GLWidget  ( QGLFormat ( QGL::DoubleBuffer    |
                                                    QGL::DepthBuffer     |
                                                    QGL::Rgba            |
                                                    QGL::AlphaChannel    |
                                                    QGL::AccumBuffer     |
                                                    QGL::StencilBuffer   |
                                                    QGL::HasOverlay      |
                                                    QGL::SampleBuffers ) ,
                                                  parent ,        p      )
{
  distance   = -40.0f              ;
  molecule   = NULL                ;
  flagBox    = false               ;
  flagAxis   = false               ;
  flagStyle  = STYLE_CPK_SPHERES   ;
  xRot       = 0                   ;
  yRot       = 0                   ;
  zRot       = 0                   ;
  background = QColor(0,0,0)       ;
  quadric    = gluNewQuadric ( )   ;
  setCursor ( Qt::OpenHandCursor ) ;
}

N::ProteinBank::PdbViewer::~PdbViewer(void)
{
  gluDeleteQuadric ( quadric ) ;
}

QSize N::ProteinBank::PdbViewer::sizeHint (void) const
{
  return QSize ( 640 , 480 ) ;
}

void N::ProteinBank::PdbViewer::setMolecule(const Molecule * molecule)
{
  this->molecule = molecule                      ;
  if (IsNull(this->molecule)) return             ;
  ////////////////////////////////////////////////
  this -> links . clear ( )                      ;
  xMin = yMin = zMin = FLT_MAX                   ;
  xMax = yMax = zMax = FLT_MIN                   ;
  ////////////////////////////////////////////////
  for (int i = 0; i < molecule->numAtoms(); i++) {
    const Atom & atom = molecule->getAtom(i)     ;
    if (atom.x() < xMin) xMin = atom . x ( )     ;
    if (atom.y() < yMin) yMin = atom . y ( )     ;
    if (atom.z() < zMin) zMin = atom . z ( )     ;
    if (atom.x() > xMax) xMax = atom . x ( )     ;
    if (atom.y() > yMax) yMax = atom . y ( )     ;
    if (atom.z() > zMax) zMax = atom . z ( )     ;
  }                                              ;
  ////////////////////////////////////////////////
  for (int i = 0; i < molecule->numLinks(); i++) {
    int id1 = molecule->getLink(i).first         ;
    int id2 = molecule->getLink(i).second        ;
    //////////////////////////////////////////////
    float distance                               ;
    distance = ::Distance                        (
               molecule->getAtom(id1).values()   ,
               molecule->getAtom(id2).values() ) ;
    //////////////////////////////////////////////
    const Atom & atom1 = molecule->getAtom(id1)  ;
    const Atom & atom2 = molecule->getAtom(id2)  ;
    int colorId1 = Atoms[atom1.code()].cpkColor  ;
    int colorId2 = Atoms[atom2.code()].cpkColor  ;
    Link link1, link2                            ;
    //////////////////////////////////////////////
    if (colorId1 == colorId2)                    {
      link1.atomId = id1                         ;
      link1.height = distance                    ;
      link1.angle  = getRotationAngle(atom1, atom2, link1.rotVec);
      this->links.push_back(link1)               ;
    } else                                       {
      link1.atomId = id1                         ;
      link1.height = distance / 2.0f             ;
      link1.angle = getRotationAngle(atom1, atom2, link1.rotVec);
      this->links.push_back(link1)               ;
      link2.atomId = id2                         ;
      link2.height = distance / 2.0f             ;
      link2.angle = getRotationAngle(atom2, atom1, link2.rotVec);
      this->links.push_back(link2)               ;
    }                                            ;
  }                                              ;
  ////////////////////////////////////////////////
  resizeGL ( width() , height () )               ;
  updateGL (                     )               ;
}

void N::ProteinBank::PdbViewer::showBox(bool enabled)
{
  flagBox = enabled ;
  updateGL ( )      ;
}

void N::ProteinBank::PdbViewer::showAxis(bool enabled)
{
  flagAxis = enabled ;
  updateGL ( )       ;
}

void N::ProteinBank::PdbViewer::initializeGL(void)
{
//  static const GLfloat lightPos[4] = { 5.0f, 5.0f, 10.0f, 1.0f } ;
//  ::glLightfv      ( GL_LIGHT0 , GL_POSITION , lightPos ) ;
  ::glEnable       ( GL_LIGHTING                        ) ;
  ::glEnable       ( GL_LIGHT0                          ) ;
  ::glEnable       ( GL_DEPTH_TEST                      ) ;
  ::glEnable       ( GL_NORMALIZE                       ) ;
  ::glLoadIdentity (                                    ) ;
  ::glTranslatef   ( 0.0f , 0.0f , distance             ) ;
}

void N::ProteinBank::PdbViewer::paintGL(void)
{
  qglClearColor  ( background                      ) ;
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) ;
  if (IsNull(molecule)) return                       ;
  glPushMatrix (                        )            ;
  glRotatef    ( xRot , 1.0 , 0.0 , 0.0 )            ;
  glRotatef    ( yRot , 0.0 , 1.0 , 0.0 )            ;
  if (flagStyle != STYLE_WIRES)                      {
    renderAtoms ( )                                  ;
    if (flagStyle != STYLE_CPK_SPHERES)              {
      renderLinks ( )                                ;
    }                                                ;
  } else                                             {
    renderWires ( )                                  ;
  }                                                  ;
  ////////////////////////////////////////////////////
  if ( flagBox  ) renderBox  ( )                     ;
  if ( flagAxis ) renderAxis ( )                     ;
  ////////////////////////////////////////////////////
  glPopMatrix ( )                                    ;
}

void N::ProteinBank::PdbViewer::resizeGL(int width,int height)
{
  double side = qMin ( width , height )     ;
  ///////////////////////////////////////////
  glViewport     ( 0 , 0 , width , height ) ;
  glMatrixMode   ( GL_PROJECTION          ) ;
  glLoadIdentity (                        ) ;
  glFrustum      ( - width  / side          ,
                     width  / side          ,
                   - height / side          ,
                     height / side          ,
                     1.0                    ,
                     2000000.0            ) ;
  glMatrixMode   ( GL_MODELVIEW           ) ;
}

float N::ProteinBank::PdbViewer::getRotationAngle (
        const AtomicPoint & coord1                ,
        const AtomicPoint & coord2                ,
        float             * vrot                  )
{
  static const float zAxis [ 3 ] = { 0.0, 0.0, 1.0 } ;
  static       float diff  [ 3 ]                     ;
  ////////////////////////////////////////////////////
  diff[0] = coord2 . x ( ) - coord1 . x ( )          ;
  diff[1] = coord2 . y ( ) - coord1 . y ( )          ;
  diff[2] = coord2 . z ( ) - coord1 . z ( )          ;
  ////////////////////////////////////////////////////
  float angle = ::vectorialAngle ( zAxis , diff )    ;
  ::vectorialProduct ( zAxis , diff , vrot )         ;
  ////////////////////////////////////////////////////
  return angle                                       ;
}

void N::ProteinBank::PdbViewer::renderAtoms(void)
{
  for (int i = 0; i < molecule->numAtoms(); i++)               {
    const Atom & atom    = molecule->getAtom(i)                ;
    int          colorId = Atoms[atom.code()].cpkColor         ;
    float        radio   = Atoms[atom.code()].vanDerWaalsRadio ;
    ////////////////////////////////////////////////////////////
    switch ( flagStyle )                                       {
      case STYLE_CPK_SPHERES                                   :
        radio /= 250.0f                                        ;
      break                                                    ;
      case STYLE_SPHERES_CYLINDERS                             :
        radio /= 1000.0f                                       ;
      break                                                    ;
      case STYLE_CYLINDERS                                     :
        radio = LINK_RADIO                                     ;
      break                                                    ;
      case STYLE_WIRES                                         :
      break                                                    ;
    }                                                          ;
    ////////////////////////////////////////////////////////////
    ::glPushMatrix (                                         ) ;
    ::glTranslatef ( atom.x() , atom.y() , atom.z()          ) ;
    ::glMaterialfv ( GL_FRONT                                  ,
                     GL_AMBIENT_AND_DIFFUSE                    ,
                     cpkTable[colorId]                       ) ;
    ::gluSphere    ( quadric , radio , 20 , 20               ) ;
    ::glPopMatrix  (                                         ) ;
    ////////////////////////////////////////////////////////////
  }                                                            ;
}

void N::ProteinBank::PdbViewer::renderLinks(void)
{
  for (int i = 0; i < links.size(); i++)                  {
    Link       & link    = links[i]                       ;
    const Atom & atom    = molecule->getAtom(link.atomId) ;
    int          colorId = Atoms[atom.code()].cpkColor    ;
    ///////////////////////////////////////////////////////
    glPushMatrix (                                      ) ;
    glTranslatef ( atom.x() , atom.y() , atom.z()       ) ;
    glRotatef    ( link . angle                           ,
                   link . rotVec [ 0 ]                    ,
                   link . rotVec [ 1 ]                    ,
                   link . rotVec [ 2 ]                  ) ;
    glMaterialfv ( GL_FRONT                               ,
                   GL_AMBIENT_AND_DIFFUSE                 ,
                   cpkTable[colorId]                    ) ;
    gluCylinder  ( quadric                                ,
                   LINK_RADIO                             ,
                   LINK_RADIO                             ,
                   link.height                            ,
                   10                                     ,
                   10                                   ) ;
    glPopMatrix  (                                      ) ;
  }                                                       ;
}

void N::ProteinBank::PdbViewer::renderWires(void)
{
  glDisable ( GL_LIGHTING )                               ;
//  glBegin   ( GL_LINES    )                               ;
  /////////////////////////////////////////////////////////
  for (int i = 0; i < links.size(); i++)                  {
    Link       & link    = links[i]                       ;
    const Atom & atom    = molecule->getAtom(link.atomId) ;
    int          colorId = Atoms[atom.code()].cpkColor    ;
    ///////////////////////////////////////////////////////
    glPushMatrix (                                      ) ;
    glTranslatef ( atom.x() , atom.y() , atom.z()       ) ;
    glRotatef    ( link.angle                             ,
                   link.rotVec[0]                         ,
                   link.rotVec[1]                         ,
                   link.rotVec[2]                       ) ;
    glBegin      ( GL_LINES                             ) ;
    glColor3fv   ( cpkTable[colorId]                    ) ;
    glVertex3d   ( 0, 0, 0                              ) ;
    glVertex3f   ( 0, 0, link.height                    ) ;
    glEnd        (                                      ) ;
    glPopMatrix  (                                      ) ;
  }                                                       ;
  /////////////////////////////////////////////////////////
//  glEnd    (             )                                ;
  glEnable ( GL_LIGHTING )                                ;
}

void N::ProteinBank::PdbViewer::renderBox(void)
{
  glDisable    ( GL_LIGHTING        ) ;
  glColor3ub   ( 255 , 255 , 255    ) ;
  glBegin      ( GL_LINES           ) ;

    glVertex3f ( xMin , yMin , zMin ) ;
    glVertex3f ( xMin , yMax , zMin ) ;
    glVertex3f ( xMin , yMax , zMin ) ;
    glVertex3f ( xMax , yMax , zMin ) ;
    glVertex3f ( xMax , yMax , zMin ) ;
    glVertex3f ( xMax , yMin , zMin ) ;
    glVertex3f ( xMax , yMin , zMin ) ;
    glVertex3f ( xMin , yMin , zMin ) ;

    glVertex3f ( xMax , yMax , zMax ) ;
    glVertex3f ( xMax , yMin , zMax ) ;
    glVertex3f ( xMax , yMin , zMax ) ;
    glVertex3f ( xMin , yMin , zMax ) ;
    glVertex3f ( xMin , yMin , zMax ) ;
    glVertex3f ( xMin , yMax , zMax ) ;
    glVertex3f ( xMin , yMax , zMax ) ;
    glVertex3f ( xMax , yMax , zMax ) ;

    glVertex3f ( xMin , yMin , zMin ) ;
    glVertex3f ( xMin , yMin , zMax ) ;
    glVertex3f ( xMin , yMax , zMin ) ;
    glVertex3f ( xMin , yMax , zMax ) ;
    glVertex3f ( xMax , yMax , zMin ) ;
    glVertex3f ( xMax , yMax , zMax ) ;
    glVertex3f ( xMax , yMin , zMin ) ;
    glVertex3f ( xMax , yMin , zMax ) ;

  glEnd        (                    ) ;
  glEnable     ( GL_LIGHTING        ) ;
}

void N::ProteinBank::PdbViewer::renderAxis(void)
{
  double xMid = (xMin + xMax) / 2   ;
  double yMid = (yMin + yMax) / 2   ;
  double zMid = (zMin + zMax) / 2   ;
  ///////////////////////////////////
  glDisable  ( GL_LIGHTING        ) ;
  glBegin    ( GL_LINES           ) ;
  ///////////////////////////////////
  glColor3ub (  255 ,    0 ,    0 ) ;
  glVertex3f ( xMin , yMid , zMid ) ;
  glVertex3f ( xMax , yMid , zMid ) ;
  ///////////////////////////////////
  glColor3ub (    0 ,  255 ,    0 ) ;
  glVertex3f ( xMid , yMin , zMid ) ;
  glVertex3f ( xMid , yMax , zMid ) ;
  ///////////////////////////////////
  glColor3ub (    0 ,    0 ,  255 ) ;
  glVertex3f ( xMid , yMid , zMin ) ;
  glVertex3f ( xMid , yMid , zMax ) ;
  ///////////////////////////////////
  glEnd      (                    ) ;
  ///////////////////////////////////
  glEnable   ( GL_LIGHTING        ) ;
}

void N::ProteinBank::PdbViewer::focusInEvent(QFocusEvent * event)
{
  if (!focusIn (event)) QGLWidget::focusInEvent (event) ;
}

void N::ProteinBank::PdbViewer::focusOutEvent(QFocusEvent * event)
{
  if (!focusOut(event)) QGLWidget::focusOutEvent(event) ;
}

void N::ProteinBank::PdbViewer::contextMenuEvent(QContextMenuEvent * event)
{
  if (Menu(event->pos())) event->accept(); else
  QGLWidget::contextMenuEvent(event);
}

void N::ProteinBank::PdbViewer::mousePressEvent(QMouseEvent * event)
{
  if (event->buttons() & Qt::LeftButton) {
    setCursor ( Qt::ClosedHandCursor )   ;
    lastPos = event -> pos ( )           ;
  }                                      ;
}

void N::ProteinBank::PdbViewer::mouseMoveEvent(QMouseEvent * event)
{
  if (event->buttons() & Qt::LeftButton) {
    int dx = event->x() - lastPos.x()    ;
    int dy = event->y() - lastPos.y()    ;
    ( xRot += dy + 360 ) %= 360          ;
    ( yRot += dx + 360 ) %= 360          ;
    lastPos = event->pos()               ;
    updateGL ( )                         ;
  }                                      ;
}

void N::ProteinBank::PdbViewer::mouseReleaseEvent(QMouseEvent * event)
{
  if (event->button() == Qt::LeftButton) {
    setCursor ( Qt::OpenHandCursor )     ;
  }                                      ;
}

void N::ProteinBank::PdbViewer::wheelEvent(QWheelEvent * event)
{
  distance -= event->delta() / 50.0f        ;
  glLoadIdentity (                        ) ;
  glTranslatef   ( 0.0f , 0.0f , distance ) ;
  updateGL       (                        ) ;
}

void N::ProteinBank::PdbViewer::setStyle(Style style)
{
  if (flagStyle != style) {
    flagStyle = style     ;
    updateGL ( )          ;
  }                       ;
}

void N::ProteinBank::PdbViewer::keyPressEvent(QKeyEvent *event)
{
}

void N::ProteinBank::PdbViewer::resizeEvent(QResizeEvent * event)
{
  QGLWidget :: resizeEvent ( event ) ;
  Relocation               (       ) ;
}

void N::ProteinBank::PdbViewer::showEvent(QShowEvent * event)
{
  QGLWidget :: showEvent ( event ) ;
  Relocation             (       ) ;
}

void N::ProteinBank::PdbViewer::Relocation(void)
{
  resizeGL ( width() , height() ) ;
}

void N::ProteinBank::PdbViewer::Export(void)
{
  QString filename = QString("%1.png").arg(windowTitle()) ;
  filename = QFileDialog::getSaveFileName                 (
               this                                       ,
               tr("Save image")                           ,
               plan->Temporary(filename)                  ,
               tr("Portable Network Graphics (*.png)")  ) ;
  if (filename.length()<=0) return                        ;
  QPixmap P = renderPixmap ( )                            ;
  QImage  I = P.toImage()                                 ;
  I . save ( filename )                                   ;
}

void N::ProteinBank::PdbViewer::startup(void)
{
  QSize ss = sizeHint ( )               ;
  resizeGL ( ss.width() , ss.height() ) ;
}

void N::ProteinBank::PdbViewer::distanceChanged(double v)
{
  distance = v ;
  updateGL ( ) ;
}

bool N::ProteinBank::PdbViewer::Menu(QPoint pos)
{ Q_UNUSED     ( pos                                ) ;
  nScopedMenu  ( mm , this                          ) ;
  QAction        * aa                                 ;
  QMenu          * ms                                 ;
  QMenu          * mb                                 ;
  QDoubleSpinBox * sb = new QDoubleSpinBox ( )        ;
  /////////////////////////////////////////////////////
  sb -> setRange    (-1000000,1000000               ) ;
  sb -> setValue    (distance                       ) ;
  sb -> setAlignment(Qt::AlignRight|Qt::AlignVCenter) ;
  /////////////////////////////////////////////////////
  mm . add     ( 101 , tr("Save image")             ) ;
  mm . addSeparator ( )                               ;
  mm . add     ( 201 , tr("Background color")       ) ;
  mm . add     ( 202 , sb                           ) ;
  mm . addSeparator ( )                               ;
  mb = mm . addMenu ( tr("Grids")                   ) ;
  mm . add     ( mb , 301 , tr("Axis") , true , flagAxis ) ;
  mm . add     ( mb , 302 , tr("Box" ) , true , flagBox  ) ;
  ms = mm . addMenu ( tr("Style")                   ) ;
  mm . add     ( ms , 401,tr("Spheres"              ),true,flagStyle == STYLE_CPK_SPHERES       ) ;
  mm . add     ( ms , 402,tr("Spheres and Cylinders"),true,flagStyle == STYLE_SPHERES_CYLINDERS ) ;
  mm . add     ( ms , 403,tr("Cylinders"            ),true,flagStyle == STYLE_CYLINDERS         ) ;
  mm . add     ( ms , 404,tr("Wires"                ),true,flagStyle == STYLE_WIRES             ) ;
  mm . setFont ( plan                               ) ;
  /////////////////////////////////////////////////////
  nConnect ( sb   , SIGNAL(valueChanged   (double))   ,
             this , SLOT  (distanceChanged(double)) ) ;
  /////////////////////////////////////////////////////
  aa = mm.exec (                                    ) ;
  if (IsNull(aa)) return true                         ;
  switch ( mm [ aa ] )                                {
    case 101                                          :
      Export ( )                                      ;
    break                                             ;
    case 201                                          :
      background = QColorDialog::getColor             (
                     background                       ,
                     this                             ,
                     tr("Background color")         ) ;
      updateGL ( )                                    ;
    break                                             ;
    case 301                                          :
      flagAxis = aa->isChecked()                      ;
      updateGL ( )                                    ;
    break                                             ;
    case 302                                          :
      flagBox  = aa->isChecked()                      ;
      updateGL ( )                                    ;
    break                                             ;
    case 401                                          :
      flagStyle  = STYLE_CPK_SPHERES                  ;
      updateGL ( )                                    ;
    break                                             ;
    case 402                                          :
      flagStyle  = STYLE_SPHERES_CYLINDERS            ;
      updateGL ( )                                    ;
    break                                             ;
    case 403                                          :
      flagStyle  = STYLE_CYLINDERS                    ;
      updateGL ( )                                    ;
    break                                             ;
    case 404                                          :
      flagStyle  = STYLE_WIRES                        ;
      updateGL ( )                                    ;
    break                                             ;
  }                                                   ;
  return true                                         ;
}

#endif
