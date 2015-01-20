import sbt._
import Keys._

object ExtendedBuild extends Build {
  lazy val versionReport = TaskKey[String]("version-report")
}
