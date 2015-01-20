
organization := "org.mbari"

name := "sbttest"

version := "1"

scalaVersion in ThisBuild := "2.11.2"

resolvers += Resolver.mavenLocal

fork := true

libraryDependencies ++= {
    val ncVersion = "4.3.23"
    Seq(
      "edu.ucar" % "netcdf" % ncVersion,
      "edu.ucar" % "grib" % ncVersion,
      "edu.ucar" % "opendap" % ncVersion,
      "edu.ucar" % "bufr" % ncVersion)
}

// CUSTOM SETTINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// set the prompt (for this build) to include the project id.
shellPrompt in ThisBuild := { state => 
  val user = System.getProperty("user.name")
  s"\n${user}@${Project.extract(state).currentRef.project}\nsbt> " 
}

// Add this setting to your project to generate a version report (See ExtendedBuild.scala too.)
// Use as 'sbt versionReport' or 'sbt version-report'
versionReport <<= (externalDependencyClasspath in Compile, streams) map {
  (cp: Seq[Attributed[File]], streams) =>
    val report = cp.map {
      attributed =>
        attributed.get(Keys.moduleID.key) match {
          case Some(moduleId) => "%40s %20s %10s %10s".format(
            moduleId.organization,
            moduleId.name,
            moduleId.revision,
            moduleId.configurations.getOrElse("")
          )
          case None           =>
            // unmanaged JAR, just
            attributed.data.getAbsolutePath
        }
    }.sortBy(a => a.trim.split("\\s+").map(_.toUpperCase).take(2).last).mkString("\n")
    streams.log.info(report)
    report
}


// For sbt-pack
packAutoSettings


// fork a new JVM for 'run' and 'test:run'
fork := true
