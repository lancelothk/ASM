package edu.cwru.cbc.ASM.commons;

import org.apache.commons.cli.*;

/**
 * Created by lancelothk on 9/18/17.
 */
public class CMDHelper {
	private String[] args;
	private String cmd_interface;
	private Options optionsToPrint;

	public CMDHelper(String[] args, String cmd_interface, Options optionsToPrint) {
		this.args = args;
		this.cmd_interface = cmd_interface;
		this.optionsToPrint = optionsToPrint;
	}

	public CommandLine build() throws ParseException {
		Options helpOption = new Options();
		helpOption.addOption(Option.builder("h").desc("Help").build());
		CommandLineParser parser = new DefaultParser();
		CommandLine helpCmd = parser.parse(helpOption, args, true);
		if (helpCmd.hasOption("h")) {
			printHelpInfo();
		}

		CommandLine cmd = null;
		try {
			cmd = parser.parse(optionsToPrint, args);
		} catch (ParseException e) {
			System.err.println(e.getMessage());
			printHelpInfo();
		}
		return cmd;
	}

	private void printHelpInfo() {
		// print the help or the version there.
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(cmd_interface, optionsToPrint);
		System.exit(Constant.EXIT_ON_SUCCESS);
	}

}
