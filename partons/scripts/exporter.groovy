//===============================================
// Analysis Script for Groovy.
// author : gavalian , date : 08/01/2018
//===============================================
import org.jlab.groot.data.*;
import org.jlab.groot.ui.*;
import org.jlab.jnp.readers.*;

String    inputFile = args[0];
Integer  dataColumn = Integer.parseInt(args[1]);

Integer xBins = 100;
Integer yBins = 100;
Integer zBins = 100;

H3F  H3 = new H3F(xBins, 0.0, 1.0, yBins, 0.0, 1.0, zBins, 0.0, 1.0);

TextFileReader reader = new TextFileReader();
reader.open(inputFile);

for(int x = 0; x < xBins; x++){
  for(int y = 0; y < yBins; y++){
    for(int z = 0; z < zBins; z++){
      reader.readNext();
      double[] array = reader.getAsDoubleArray();
      //double value = array[dataColumn];
      double value = reader.getAsDouble(dataColumn);
      System.out.println(value);
      H3.setBinContent(x,y,z,value);
    }
  }
}

H3.export("data.nrrd");
