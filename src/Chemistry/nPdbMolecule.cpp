#include <chemistry.h>

float parseFloat(const QString & str,int pos,int len)
{
  float value = -1.0                    ;
  bool  ok                              ;
  if (pos + len <= str.length())        {
    QStringRef ref(&str, pos, len)      ;
    value = ref.toString().toFloat(&ok) ;
  }                                     ;
  return value                          ;
}

int parseInt(const QString & str,int pos,int len)
{
  bool ok                               ;
  int  value = -1                       ;
  if (pos <= str.length())              {
    len = qMin(len, str.length() - pos) ;
    QStringRef ref ( &str , pos , len ) ;
    value = ref.toString().toInt(&ok)   ;
    if (!ok) value = -1                 ;
  }                                     ;
  return value                          ;
}

//////////////////////////////////////////////////////////////////////////

N::ProteinBank::Molecule:: Molecule(QObject * parent)
                         : QObject (          parent)
{
}

N::ProteinBank::Molecule::~Molecule(void)
{
}

void N::ProteinBank::Molecule::load(const QString &filename)
{
  QFile file ( filename )                 ;
  if (!file.open(QFile::ReadOnly)) return ;
  QTextStream input ( &file )             ;
  parse ( input )                         ;
}

void N::ProteinBank::Molecule::load(const QByteArray & body)
{
  QTextStream input ( body  ) ;
  parse             ( input ) ;
}

void N::ProteinBank::Molecule::parse(QTextStream & input)
{
  QString line                                                ;
  m_Atoms . clear ( )                                         ;
  m_Links . clear ( )                                         ;
  while (line = input.readLine(), !line.isNull())             {
    if (line.startsWith("ATOM") || line.startsWith("HETATM")) {
      Atom    atom                                            ;
      QString codeStr = QStringRef(&line, 12, 4).toString()   ;
      atom . setX    ( parseFloat(line, 30, 8) )              ;
      atom . setY    ( parseFloat(line, 38, 8) )              ;
      atom . setZ    ( parseFloat(line, 46, 8) )              ;
      atom . setCode ( codeStr.trimmed()       )              ;
      m_Atoms . push_back ( atom )                            ;
    } else
    if (line.startsWith("CONECT"))                            {
      int id1 = parseInt ( line ,  6 , 5 )                    ;
      int id2 = parseInt ( line , 11 , 5 )                    ;
      int id3 = parseInt ( line , 16 , 5 )                    ;
      int id4 = parseInt ( line , 21 , 5 )                    ;
      int id5 = parseInt ( line , 26 , 5 )                    ;
      if (id1 != -1)                                          {
        if (id2 != -1) addLink ( id1 - 1 , id2 - 1 )          ;
        if (id3 != -1) addLink ( id1 - 1 , id3 - 1 )          ;
        if (id4 != -1) addLink ( id1 - 1 , id4 - 1 )          ;
        if (id5 != -1) addLink ( id1 - 1 , id5 - 1 )          ;
      }                                                       ;
    }                                                         ;
  }                                                           ;
}

int N::ProteinBank::Molecule::numAtoms(void) const
{
  return m_Atoms . size ( ) ;
}

int N::ProteinBank::Molecule::numLinks(void) const
{
  return m_Links . size ( ) ;
}

const N::ProteinBank::Atom & N::ProteinBank::Molecule::getAtom(int id) const
{
  return m_Atoms [ id ] ;
}

const N::ProteinBank::Molecule::Link & N::ProteinBank::Molecule::getLink(int id) const
{
  return m_Links[id];
}

void N::ProteinBank::Molecule::addLink(int id1,int id2)
{
  if ( id1 > id2 ) qSwap ( id1 , id2 )                             ;
  Link newEnlace ( id1 , id2 )                                     ;
  //////////////////////////////////////////////////////////////////
  QVector<Link>::const_iterator first,last                         ;
  first = m_Links . begin ( )                                      ;
  last  = m_Links . end   ( )                                      ;
  //////////////////////////////////////////////////////////////////
  int numAtomos = this -> numAtoms ( )                             ;
  if ( 0 <= id1 && id1 < numAtomos && 0 <= id2 && id2 < numAtomos) {
    if (qFind(first, last, newEnlace) == last)                     {
      m_Links . push_back ( Link ( id1 , id2 ) )                   ;
    }                                                              ;
  }                                                                ;
}
