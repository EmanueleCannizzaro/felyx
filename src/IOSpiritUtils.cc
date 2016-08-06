  #include "IOSpiritUtils.h"


  unsigned long fe_base::ReadFileToMemory(const char* filepath, char* &buffer)
  {
    unsigned long size = 0;

    std::ifstream infile(filepath);

    if (!infile.good())
    {
      std::string error("ERROR: In ReadFileToMemory: Could not open the file ");
      error += filepath;
      FELYX_RUNTIME_THROW( error.c_str() );
    }

    // get size of file
    infile.seekg(0,std::ifstream::end);
    size=infile.tellg();
    infile.seekg(0);

    // allocate memory for file content
    buffer = new char [size];

    // read content of infile
    infile.read (buffer,size);

    infile.close();

    return size;
  }
