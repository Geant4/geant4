//GGE (GEANT4 Geometry Editor)
//GGE Item
//Tetsuya Yamada

abstract class GGEItem implements java.io.Serializable {
	public String name;
	public boolean referenced = false;

	public int compareName(GGEItem gi) {
	  return name.compareTo(gi.name);
	}

}