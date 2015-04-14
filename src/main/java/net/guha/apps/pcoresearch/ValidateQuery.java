package net.guha.apps.pcoresearch;

import org.xml.sax.SAXException;

import javax.xml.XMLConstants;
import javax.xml.transform.Source;
import javax.xml.transform.stream.StreamSource;
import javax.xml.validation.Schema;
import javax.xml.validation.SchemaFactory;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;

/**
 * Validate an XML pharmacophore definition file.
 * <p/>
 * See http://stackoverflow.com/questions/1541253/how-to-validate-an-xml-document-using-a-relax-ng-schema-and-jaxp for
 * the description of how to use JING to validate Relax-NG compact schemata. Seems to work fine on JDK 1.6
 *
 * @author Rajarshi Guha
 */
public class ValidateQuery {
    javax.xml.validation.Validator validator;

    /**
     * Create validator using built in RelaxNG schema definition.
     *
     * @throws SAXException
     */
    public ValidateQuery() throws SAXException {
        System.setProperty(SchemaFactory.class.getName() + ":" + XMLConstants.RELAXNG_NS_URI, "com.thaiopensource.relaxng.jaxp.XMLSyntaxSchemaFactory");
        SchemaFactory factory = SchemaFactory.newInstance(XMLConstants.RELAXNG_NS_URI);
        InputStream is = getClass().getResourceAsStream("/pharmacophore.rng");
        Schema schema = factory.newSchema(new StreamSource(is));
        validator = schema.newValidator();
    }

    /**
     * Create validator using supplied RelaxNG schema definition.
     *
     * @param schemaFileName The file name of the schema
     * @throws SAXException
     */
    public ValidateQuery(String schemaFileName) throws SAXException {
        System.setProperty(SchemaFactory.class.getName() + ":" + XMLConstants.RELAXNG_NS_URI, "com.thaiopensource.relaxng.jaxp.XMLSyntaxSchemaFactory");
        SchemaFactory factory = SchemaFactory.newInstance(XMLConstants.RELAXNG_NS_URI);
        Schema schema = factory.newSchema(new File(schemaFileName));
        validator = schema.newValidator();
    }


    public boolean validate(String fileName) throws IOException, SAXException {
        Source source = new StreamSource(fileName);
        try {
            validator.validate(source);
            return true;
        } catch (SAXException ex) {
            throw new SAXException(ex.getMessage());
        }
    }

    public static void main(String[] args) throws IOException, SAXException {
        ValidateQuery vq = new ValidateQuery("/Users/guhar/src/cdkpcore/pharmacophore.rng");
        vq.validate("/Users/guhar/src/cdkpcore/data/simple.xml");
    }


}
