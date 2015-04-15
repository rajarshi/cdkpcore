CDKPsearch is a command line to tool to perform simple pharmacophore matching based on a pharmacophore definition. To build the tool
you will need [Maven](https://maven.apache.org/). The executable JAR file can be built using the following command
```
mvn clean package
```
The resultant jar can be found in the ```target/``` folder

Usage
-----
The tool performs the search  on a collection of molecules stored in SD format. The program accepts an SD file with single (i.e., one conformer) structures or multi-conformer structures. In the latter case, conformers are detected based on titles. Thus all conformers for a given molecule should be located in sequence and should have the same title. The program will write out the structures that match the query to the file ```hits.sdf``` and also provide a summary report in ```report.txt```. You can run the program as
```
java -jar CDKPsearch-1.3.0.jar --sdfile targets.sdf --query query.xml -c 
```

Pharmacophore Definition Format
-------------------------------

The format for the pharmacophore query file is XML with a [RelaxNG](http://relaxng.org/) schema available [here](https://github.com/rajarshi/cdkpcore/blob/master/src/main/resources/pharmacophore.rng). This can be used to validate any pharmacophore definition file using tools such as [jing](http://www.thaiopensource.com/relaxng/jing.html). An example of a query file containng multiple pharmacophore queries is given below:
```xml
<?xml version="1.0"?>
<pharmacophoreContainer version="1.0">
    <!--
     some gobal pharmacphore group definitions, usable by 
        all pharmacophores defined in this file 
    -->
    <group id="Ha" description="H-bond acceptor">
        [!$([#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]
    </group>
    <group id="Hd" description="H-bond donor">[!$([#6,H0,-,-2,-3])]</group>
    <group id="Aromatic">[*;a]</group>
    <group id="Acid" description="Acidic groups, defined as a proton donor">[!H0;F,Cl,Br,I,N+,$([OH]-*=[!#6]),+]</group>
    <pharmacophore description="An imaginary pharmacophore definition" name="Pcore Def 1">
        <!--
         in this case we only want to consider carboxylic acid groups
        	    so we provide a local definition
        	
        -->
        <group id="CarbAcid">[CX3](=O)[OX2H1]</group>
        <!--  An exact distance constraint  -->
        <distanceConstraint lower="1.4" units="A">
            <groupRef id="Hd"/>
            <groupRef id="CarbAcid"/>
        </distanceConstraint>
        <!--  A distance range constraint  -->
        <distanceConstraint lower="3.5" upper="4.8" units="A">
            <groupRef id="Ha"/>
            <groupRef id="Acid"/>
        </distanceConstraint>
        <!--   An angle range constraint  -->
        <angleConstraint lower="46" upper="47.5" units="degrees">
            <groupRef id="Ha"/>
            <groupRef id="Acid"/>
            <groupRef id="Ha"/>
        </angleConstraint>
    </pharmacophore>
    <pharmacophore description="A definition for the D1 receptor" name="D1">
        <group id="Hydroxyl">[OX2H]</group>
        <group id="BasicAmine">[NX3;H2,H1;!$(NC=O)]</group>
        <distanceConstraint lower="2.7" upper="2.9" units="A">
            <groupRef id="Aromatic"/>
            <groupRef id="Hydroxyl"/>
        </distanceConstraint>
        <distanceConstraint lower="4.2" upper="4.8" units="A">
            <groupRef id="Aromatic"/>
            <groupRef id="BasicAmine"/>
        </distanceConstraint>
        <distanceConstraint lower="6.8" upper="8.3" units="A">
            <groupRef id="Hydroxyl"/>
            <groupRef id="BasicAmine"/>
        </distanceConstraint>
    </pharmacophore>
</pharmacophoreContainer>
```
Note that the CDKPsearch application assumes that a query file will have a single query and if multiple queries are present will complain. This will be updated in the future. The schema supports both distance and angle constraint, both are currently supported by the CDK. Dihedral constraints are on the way.
