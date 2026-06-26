// PythiaStdlib.cc is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for Pythia utilities.

#include "Pythia8/PythiaStdlib.h"

// Directory creation for POSIX.
#include <sys/stat.h>

namespace Pythia8 {

//==========================================================================

// Convert string to lowercase for case-insensitive comparisons.
// By default remove any initial and trailing blanks or escape characters.

string toLower(const string& name, bool trim) {

  // Copy string without initial and trailing blanks or escape characters.
  string temp = name;
  if (trim) temp = trimString(name);

  // Convert to lowercase letter by letter.
  for (int i = 0; i < int(temp.length()); ++i) temp[i] = tolower(temp[i]);
  return temp;

}

//==========================================================================

// Return the given string with any initial and trailing blanks or
// escape characters removed.

string trimString(const string& name) {

  // Copy string without initial and trailing blanks or escape characters.
  if (name.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return "";
  int firstChar = name.find_first_not_of(" \n\t\v\b\r\f\a");
  int lastChar  = name.find_last_not_of(" \n\t\v\b\r\f\a");
  return name.substr( firstChar, lastChar + 1 - firstChar);

}

//==========================================================================

// Convert a double to string using reasonable formatting.

string toString(double val) {

  stringstream ssval;
  bool eform = false;
  if ( val == 0. ) ssval << fixed << setprecision(1);
  else if ( abs(val) < 0.001 ) ssval << scientific << setprecision(4);
  else if ( abs(val) < 0.1 ) ssval << fixed << setprecision(7);
  else if ( abs(val) < 1000. ) ssval << fixed << setprecision(5);
  else if ( abs(val) < 1000000. ) ssval << fixed << setprecision(3);
  else {eform = true; ssval << scientific << setprecision(4);}
  ssval << val;
  string sval = ssval.str();
  if (!eform) sval.erase(sval.find_last_not_of('0') + 1);
  return sval;

}

//==========================================================================

// Split a string by a delimiter.

vector<string> splitString(string val, string delim) {

  vector<string> vectorVal;
  size_t         stringPos(0);
  while (stringPos != string::npos) {
    stringPos = val.find(delim);
    vectorVal.push_back(val.substr(0, stringPos));
    val = val.substr(stringPos + 1);
  }
  return vectorVal;

}

//==========================================================================

// Allow several alternative inputs for true/false.

bool boolString(string tag) {

  string tagLow = toLower(tag);
  return ( tagLow == "true" || tagLow == "1" || tagLow == "on"
  || tagLow == "yes" || tagLow == "ok" );

}

//==========================================================================

// Extract XML value string following XML attribute.

string attributeValue(string line, string attribute) {

  if (line.find(attribute) == string::npos) return "";
  int iBegAttri = line.find(attribute);
  int iBegQuote = line.find("\"", iBegAttri + 1);
  int iEndQuote = line.find("\"", iBegQuote + 1);
  return line.substr(iBegQuote + 1, iEndQuote - iBegQuote - 1);

}

//==========================================================================

// Extract XML bool value following XML attribute.

bool boolAttributeValue(string line, string attribute) {

  string valString = attributeValue(line, attribute);
  if (valString == "") return false;
  return boolString(valString);

}

//==========================================================================

// Extract XML int value following XML attribute.

int intAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0;
  istringstream valStream(valString);
  int intVal;
  valStream >> intVal;
  return intVal;

}

//==========================================================================

// Extract XML double value following XML attribute.

double doubleAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0.;
  istringstream valStream(valString);
  double doubleVal;
  valStream >> doubleVal;
  return doubleVal;

}

//==========================================================================

// Extract XML bool vector value following XML attribute.

vector<bool> boolVectorAttributeValue(string line,
  string attribute) {
  string valString = attributeValue(line, attribute);
  size_t openBrace  = valString.find_first_of("{");
  size_t closeBrace = valString.find_last_of("}");
  if (openBrace != string::npos)
    valString = valString.substr(openBrace + 1, closeBrace - openBrace - 1);
  if (valString == "") return vector<bool>();
  vector<bool> vectorVal;
  size_t       stringPos(0);
  while (stringPos != string::npos) {
    stringPos = valString.find(",");
    istringstream  valStream(valString.substr(0, stringPos));
    valString = valString.substr(stringPos + 1);
    vectorVal.push_back(boolString(valStream.str()));
  }
  return vectorVal;

}

//==========================================================================

// Extract XML int vector value following XML attribute.

vector<int> intVectorAttributeValue(string line,
  string attribute) {
  string valString = attributeValue(line, attribute);
  size_t openBrace  = valString.find_first_of("{");
  size_t closeBrace = valString.find_last_of("}");
  if (openBrace != string::npos)
    valString = valString.substr(openBrace + 1, closeBrace - openBrace - 1);
  if (valString == "") return vector<int>();
  int         intVal;
  vector<int> vectorVal;
  size_t      stringPos(0);
  while (stringPos != string::npos) {
    stringPos = valString.find(",");
    istringstream  valStream(valString.substr(0, stringPos));
    valString = valString.substr(stringPos + 1);
    valStream >> intVal;
    vectorVal.push_back(intVal);
  }
  return vectorVal;

}

//==========================================================================

// Extract XML double vector value following XML attribute.

vector<double> doubleVectorAttributeValue(string line,
  string attribute) {
  string valString = attributeValue(line, attribute);
  size_t openBrace  = valString.find_first_of("{");
  size_t closeBrace = valString.find_last_of("}");
  if (openBrace != string::npos)
    valString = valString.substr(openBrace + 1, closeBrace - openBrace - 1);
  if (valString == "") return vector<double>();
  double         doubleVal;
  vector<double> vectorVal;
  size_t         stringPos(0);
  while (stringPos != string::npos) {
    stringPos = valString.find(",");
    istringstream  valStream(valString.substr(0, stringPos));
    valString = valString.substr(stringPos + 1);
    valStream >> doubleVal;
    vectorVal.push_back(doubleVal);
  }
  return vectorVal;

}

//==========================================================================

// Extract XML string vector value following XML attribute.

vector<string> stringVectorAttributeValue(string line,
  string attribute) {
  string valString = attributeValue(line, attribute);
  size_t openBrace  = valString.find_first_of("{");
  size_t closeBrace = valString.find_last_of("}");
  if (openBrace != string::npos)
    valString = valString.substr(openBrace + 1, closeBrace - openBrace - 1);
  if (valString == "") return vector<string>();
  string         stringVal;
  vector<string> vectorVal;
  size_t         stringPos(0);
  while (stringPos != string::npos) {
    stringPos = valString.find(",");
    if (stringPos != string::npos) {
      vectorVal.push_back(valString.substr(0, stringPos));
      valString = valString.substr(stringPos + 1);
    } else vectorVal.push_back(valString);
  }
  return vectorVal;

}

//==========================================================================

// Complete an XML tag.

void completeTag(istream& stream, string& line) {
  while (line.find(">") == string::npos) {
    string addLine;
    if (!getline(stream, addLine)) break;
    line += " " + addLine;
  }
}

//==========================================================================

} // end namespace Pythia8
