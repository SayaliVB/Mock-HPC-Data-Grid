package gash.data.nyc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import com.fasterxml.jackson.databind.MappingIterator;
import com.fasterxml.jackson.dataformat.csv.CsvMapper;
import com.fasterxml.jackson.dataformat.csv.CsvParser;

/**
 * Example data set from NYC data (https://data.cityofnewyork.us). The goal of
 * this project is to show how to parse the NYC data set and transform (ETL) the
 * data into JSON format. Note you are not required to use JSON, or this code.
 * 
 * @TODO Your challenge is to implement all code required for the transformation
 *       process to cleanly remove improper records and validate all fields.
 *       This is in preparation for the lightning challenge.
 * 
 *       <pre>
 * Test path:
 * ~/data/
 * 
 * Data Files:
 * Parking_Violations_Issued_-_Fiscal_Year_[2014-2022].csv
 * 
 * Data Header (first line):
 * Summons Number,Plate ID,Registration State,Plate Type,
 * Issue Date,Violation Code,Vehicle Body Type,Vehicle Make,
 * Issuing Agency,Street Code1,Street Code2,Street Code3,
 * Vehicle Expiration Date,Violation Location,Violation Precinct,
 * Issuer Precinct,Issuer Code,Issuer Command,Issuer Squad,
 * Violation Time,Time First Observed,Violation County,
 * Violation In Front Of Or Opposite,House Number,Street Name,
 * Intersecting Street,Date First Observed,Law Section,Sub Division,
 * Violation Legal Code,Days Parking In Effect,From Hours In Effect,
 * To Hours In Effect,Vehicle Color,Unregistered Vehicle?,Vehicle Year,
 * Meter Number,Feet From Curb,Violation Post Code,Violation Description,
 * No Standing or Stopping Violation,Hydrant Violation,Double Parking Violation
 * 
 * 
 * Codes Header (first line) File: ParkingViolationCodes_January2020.csv
 * VIOLATION CODE,VIOLATION DESCRIPTION,"Manhattan  96th St. & below
 * (Fine Amount $)","All Other Areas
 * (Fine Amount $)"

 * Code Data:
 * 46,DOUBLE PARKING,115,115
 * 47,DOUBLE PARKING-MIDTOWN COMML,115,115
 * 48,BIKE LANE,115,115
 *
 *       </pre>
 * 
 *       copyright 2023, gash
 *
 *       Gash licenses this file to you under the Apache License, version 2.0
 *       (the "License"); you may not use this file except in compliance with
 *       the License. You may obtain a copy of the License at:
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *       Unless required by applicable law or agreed to in writing, software
 *       distributed under the License is distributed on an "AS IS" BASIS,
 *       WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 *       implied. See the License for the specific language governing
 *       permissions and limitations under the License.
 * 
 * @author gash
 *
 */
public class NYCParkingParser {

	public static final String basename = "NYC_Parking_Violations_Issued";

	/** json field names */
	private static final String[] names = { "SummonsNumber", "PlateID", "RegistrationState", "PlateType", "IssueDate",
			"ViolationCode", "VehicleBodyType", "VehicleMake", "IssuingAgency", "StreetCode1", "StreetCode2",
			"StreetCode3", "VehicleExpirationDate", "ViolationLocation", "ViolationPrecinct", "IssuerPrecinct",
			"IssuerCode", "IssuerCommand", "IssuerSquad", "ViolationTime", "TimeFirstObserved", "ViolationCounty",
			"ViolationInFrontOfOrOpposite", "HouseNumber", "StreetName", "IntersectingStreet", "DateFirstObserved",
			"LawSection", "SubDivision", "ViolationLegalCode", "DaysParkingInEffect", "FromHoursInEffect",
			"ToHoursInEffect", "VehicleColor", "UnregisteredVehicle?", "VehicleYear", "MeterNumber", "FeetFromCurb",
			"ViolationPostCode", "ViolationDescription", "NoStandingorStoppingViolation", "HydrantViolation",
			"DoubleParkingViolation" };

	/**
	 * generate filename
	 * 
	 * @param yr
	 *            year (YYYY)
	 * @param mo
	 *            month (MM, 1-12)
	 * @param day
	 *            day (DD, 01-31)
	 * @param code
	 *            (code, 01-99)
	 * @return
	 */
	public static String filename(int yr, int mo, int day, int code) {
		var r = String.format("%d/%02d/%02d/%s_%02d.csv", yr, mo, day, NYCParkingParser.basename, code);
		return r;
	}

	@SuppressWarnings("unchecked")
	public static void processRawToJson(File rawFile, File outDirectory) throws Exception {
		if (!rawFile.exists())
			throw new IOException("Raw file does not exist or not accessible.");
		if (!outDirectory.exists()) {
			if (!outDirectory.mkdirs())
				throw new IOException("Output directory not accessible.");
		}

		var jsons = new HashMap<String, JSONArray>();

		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(rawFile));

			// header
			var raw = br.readLine();
			if (raw == null)
				throw new IOException("File is empty");

			var header = new JSONObject();
			header.put("description", "NYC Parking Violation Data");
			header.put("origin", rawFile.getName());
			header.put("ref", "https://data.cityofnewyork.us");
			header.put("created", System.currentTimeMillis());

			// header
			while (raw != null) {
				raw = br.readLine();
				if (raw == null)
					continue;

				// parts - does not appear to have embedded commas in fields.
				// arg -1, want empty fields! length() -1 for ^M
				var parts = raw.substring(0, raw.length() - 1).split(",", -1);

				// TODO this is to show you what is being read. Not practical
				// outside of testing!
				System.out.println("parts[0]: " + parts[0]);
				System.out.println("parts[1]: " + parts[1]);
				System.out.println("parts[2]: " + parts[2]);
				System.out.println("parts[3]: " + parts[3]);
				System.out.println("parts[4]: " + parts[4]);
				System.out.println("parts[5]: " + parts[5]);
				System.out.println("parts[6]: " + parts[6]);

				// TODO What do you want to use to organize your data?
				//
				// organize data into a binary tree
				var subfilename = String.format("%s/%s_%s.json", dateStr(parts[4]), basename, parts[5]);
				var doc = jsons.get(subfilename);
				if (doc == null) {
					doc = new JSONArray();
					jsons.put(subfilename, doc);
				}

				var entry = new JSONObject();
				for (int i = 0; i < parts.length; i++) {
					entry.put(names[i], parts[i].trim());
				}
				doc.add(entry);

			}

			// at the end of reading the raw file, write out parsed data. Note
			// some files may already exist
			for (var name : jsons.keySet()) {
				var data = jsons.get(name);
				var outF = new File(outDirectory, name);
				outF.getParentFile().mkdirs();
				var out = new BufferedWriter(new FileWriter(outF));
				if (!outF.exists()) {
					data.add(0, header);
				}

				// TODO this could be too big of a string, maybe iterate and
				// write one row at a time
				out.write(data.toJSONString());
				out.flush();
				out.close();
			}
		} finally {
			if (br != null)
				br.close();
		}
	}

	@SuppressWarnings("unchecked")
	public static void processJacksonCSVToJson(File rawFile, File outDirectory) throws Exception {
		if (!rawFile.exists())
			throw new IOException("Raw file does not exist or not accessible.");
		if (!outDirectory.exists()) {
			if (!outDirectory.mkdirs())
				throw new IOException("Output directory not accessible.");
		}

		var jsons = new HashMap<String, JSONArray>();

		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(rawFile));

			// header
			var raw = br.readLine();
			if (raw == null)
				throw new IOException("File is empty");

			var header = new JSONObject();
			header.put("description", "NYC Parking Violation Data");
			header.put("origin", rawFile.getName());
			header.put("source", rawFile.getName());
			header.put("ref", "https://data.cityofnewyork.us");

			int maxRows = 250000;
			int count = 0;

			CsvMapper mapper = new CsvMapper();
			MappingIterator<List<String>> iter = mapper.readerForListOf(String.class)
					.with(CsvParser.Feature.WRAP_AS_ARRAY).readValues(br);
			while (iter.hasNext()) {
				var row = iter.next();

				var debug = false;
				if (debug) {
					System.out.println("parts[0]: " + row.get(0));
					System.out.println("parts[1]: " + row.get(1));
					System.out.println("parts[2]: " + row.get(2));
					System.out.println("parts[3]: " + row.get(3));
					System.out.println("parts[4]: " + row.get(4));
					System.out.println("parts[5]: " + row.get(5));
					System.out.println("parts[6]: " + row.get(6));
				}

				// the key becomes the filename
				var yr = dateStr(row.get(4));
				var yrfn = yr.replace('/', '_');
				var key = String.format("%s/%s_%s_%s.json", yr, yrfn, basename, row.get(5));
				var doc = jsons.get(key);
				if (doc == null) {
					doc = new JSONArray();
					jsons.put(key, doc);
				}

				var entry = new JSONObject();
				for (int i = 0; i < row.size(); i++) {
					entry.put(names[i], row.get(i));
				}
				doc.add(entry);

				// limit memory consumption by periodically flushing to file
				count++;
				if (count > maxRows) {
					flushAndCloseLists(outDirectory, jsons, header, false);
					count = 0;
				}

			}

			flushAndCloseLists(outDirectory, jsons, header, true);
		} finally {
			if (br != null)
				br.close();
		}
	}

	private static void flushAndCloseLists(File outDirectory, HashMap<String, JSONArray> jsons, JSONObject header,
			boolean done) throws Exception {
		BufferedWriter out = null;
		try {
			// at the end of reading the raw file, write out parsed data. Note
			// some files may already exist
			for (var name : jsons.keySet()) {
				var data = jsons.get(name);
				var outF = new File(outDirectory, name);

				if (!outF.getParentFile().exists())
					outF.getParentFile().mkdirs();

				// open in append mode
				var isNew = !outF.exists();
				out = new BufferedWriter(new FileWriter(outF, true));
				if (isNew) {
					out.write("{\n");
					out.write("\"header\":" + header.toJSONString() + ",\n");
					out.write("\"data\": [\n");
					out.flush();
				}

				// TODO this could be too big of a string, maybe iterate and
				// write one row at a time
				for (var entry : data) {
					out.write(((JSONObject) entry).toJSONString() + ",\n");
				}

				if (data.size() > 0)
					System.out.println("writing " + data.size() + " rows to " + outF.getAbsolutePath());

				// close the json document and array
				if (done) {
					out.write("]\n}");
				}

				out.flush();
				out.close();
				data.clear();
			}
		} finally {
			if (out != null)
				out.close();
		}
	}

	/**
	 * use the date to create a subdirectory.
	 * @param rawDate as mo/dy/year
	 * @return
	 * @throws IllegalArgumentException
	 */
	private static String dateStr(String rawDate) throws IllegalArgumentException {
		// TODO what could go wrong?
		var parts = rawDate.split("/");
		return String.format("%s/%s/%s", parts[3], parts[0], parts[1]);
	}

	/**
	 * Example main/app to run from the command line.
	 * 
	 * Usage: java [class] source-file out-directory
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		if (args.length != 2)
			throw new RuntimeException("Usage: main <in file> <out dir>");

		var src = new File(args[0]);
		var outdir = new File(args[1]);

		if (!outdir.exists() && !outdir.mkdirs())
			throw new RuntimeException("Cannot create " + outdir.getAbsolutePath());

		if (!src.exists())
			throw new RuntimeException("Input file does not exist");

		NYCParkingParser.processJacksonCSVToJson(src, outdir);

	}
}
