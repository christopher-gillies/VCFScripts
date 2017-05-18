
############
# Classes for writing out makefile
############

class MakeEntry:
	
	def __init__(self, target, commands, dependencies, comment = ""):
		#target is the output file
		#commands should be an array of shell commands
		#dependencies should be MakeEntries dependent on this
		self.target = target
		self.commands = commands
		self.dependencies = []
		self.comment = comment
		for dep in dependencies:
			self.dependencies.append(  dep.target )
			
		self.commands.append("touch {0}".format(self.target))
		
	def __str__(self):
		commands_adjust = []
		for com in self.commands:
			commands_adjust.append("\t{0}".format(com))
		return "#{comment}\n{tar}:\t{dep}\n{com}".format(comment=self.comment, tar=self.target, dep="\t".join(self.dependencies), com="\n".join(commands_adjust))


class MakeFile:
	def __init__(self, target, entries):
		self.entries = entries
		self.dependencies = []
		self.target = target
		self.clean_entry = None
		for dep in entries:
			self.dependencies.append(  dep.target )
	
	def set_clean_entry(self, clean_entry):
		self.clean_entry = clean_entry
		
	def __str__(self):
		root = "all:\t{dep}\n\ttouch {target}".format(dep="\t".join(self.dependencies), target = self.target)
		entries_str = []
		for entry in self.entries:
			entries_str.append( str(entry) )
		
		if self.clean_entry is not None:
			entries_str.append( str(self.clean_entry ) )
		
		return "{0}\n\n{1}".format(root, "\n\n".join(entries_str))


#mk_test = MakeEntry("11111", ["print a","print b"], [])	
#mk_test_a = MakeEntry("2222", ["print c","print d"], [mk_test])
#mk_file = MakeFile("all.ok", [ mk_test,mk_test_a ])
#print mk_file
