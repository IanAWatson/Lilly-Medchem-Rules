class Bindir
  include Enumerable

  def initialize(home)
    @uname = `uname`.chomp!
    @version = `uname -r`.chomp!
    @bindir = Array.new

    set_home(home)
  end

  def add_dir(dir)

    return 0 unless FileTest.directory?(dir)

    @bindir.unshift(dir)

    return 1
  end

  def set_home home

    raise "Invalid binary directory '#{home}'" unless FileTest.directory?(home)

    rc = 0

    rc += add_dir("#{home}/bin")

    rc += add_dir("#{home}/bin/#{@uname}/")

    rc += add_dir("#{home}/bin/#{@uname}-#{@version}/")

    $stderr.print "No '#{@uname}' or '#{@uname}-#{@version}' directories found under '#{home}'\n" unless rc > 0
  end

  def push dir
    raise "Invalid binary directory '#{dir}'" unless FileTest.directory?(dir)
    @bindir.push(dir)
  end

  def report_not_found prog
    $stderr.print "Cannot find '#{prog}' in any of\n"
    @bindir.each do |dir|
      $stderr.print "#{dir}\n"
    end

    return prog
  end

  def find_executable prog
    @bindir.each { |dir|
      fname = File.join(dir, prog)

#     $stderr.print "Checking '#{fname}'\n"

      return fname if File.executable_real? fname
    }

    return report_not_found(prog)
  end

# Kind of awkward how this should be implemented. Just to
# simplify, look for the shell wrapper first, then the executable

  def find_shell_wrapper_or_executable prog
    @bindir.each do |dir|
      fname = File.join(dir, "#{prog}.sh")

#     $stderr.print "Checking '#{fname}'\n"

      return fname if File.executable_real? fname
    end

    return find_executable(prog)
  end

  def each
    @bindir.each { |dir|
      yield dir
    }
  end
end
