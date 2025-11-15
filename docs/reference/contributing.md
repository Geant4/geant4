# Contributing to Geant4

## Overview

Contributions to Geant4 are welcome through GitHub Pull Requests. This guide outlines the process and requirements for contributing.

## Contribution Guidelines

From CONTRIBUTING.rst:3-11:

### Pull Request Requirements

- **Single Topic**: A single PR should refer to a single particular topic (fix or suggested enhancement)
- **No Uncorrelated Changes**: PRs including uncorrelated changes in different areas or modules will not be considered
- **Review Process**: PRs will be examined by a responsible party and potentially extracted for inclusion in the development repository

### Pull Request Process

1. Create a PR on GitHub with your changes
2. Ensure your PR focuses on a single topic
3. The PR will be reviewed by the Geant4 Collaboration
4. If accepted, it will be included in a future release or patch

## Bug Reporting

For bug fixes, it is recommended to:

1. Report the issue through the official [Bugzilla problem reporting system](https://bugzilla-geant4.kek.jp)
2. Link back to your Pull Request if you've created one
3. Provide detailed information about the bug and how to reproduce it

From CONTRIBUTING.rst:13-16:
```
For simple bug-fixes, it is still recommended to use the official Bugzilla
problem reporting system to report the issue, linking back to the Pull
Request if created
```

## Feature Requests & Questions

- **Questions**: Post on the [User Forum](https://geant4-forum.web.cern.ch)
- **Feature Requests**: Submit through the User Forum
- **Documentation**: Review [full documentation](https://cern.ch/geant4/support/user_documentation)

From CONTRIBUTING.rst:18-19:
```
Questions on Geant4 itself or feature requests should be made through our
User Forum, and full documentation on use of the toolkit is online
```

## Code Quality Standards

When contributing code, ensure:

### Code Style

Geant4 includes code formatting configurations:
- **Clang Format**: `.clang-format` in the root directory
- **Clang Tidy**: `.clang-tidy` in the root directory

Use these tools to maintain consistent code style:

```bash
# Format your code
clang-format -i your_file.cc

# Check with clang-tidy
clang-tidy your_file.cc
```

### Best Practices

- Follow C++ best practices
- Write clear, self-documenting code
- Include comments for complex logic
- Maintain backward compatibility when possible
- Add appropriate unit tests if applicable

## Development Workflow

### 1. Fork and Clone

```bash
git clone https://github.com/YOUR_USERNAME/geant4.git
cd geant4
```

### 2. Create a Branch

```bash
git checkout -b fix/my-bug-fix
# or
git checkout -b feature/my-enhancement
```

### 3. Make Changes

- Follow the code style guidelines
- Keep changes focused on the single topic
- Test your changes thoroughly

### 4. Commit

```bash
git add .
git commit -m "Brief description of changes"
```

### 5. Push and Create PR

```bash
git push origin fix/my-bug-fix
```

Then create a Pull Request on GitHub.

## Getting Help

If you need assistance with contributions:

- **Forum**: [https://geant4-forum.web.cern.ch](https://geant4-forum.web.cern.ch)
- **Documentation**: [https://cern.ch/geant4/support/user_documentation](https://cern.ch/geant4/support/user_documentation)
- **Bug Tracker**: [https://bugzilla-geant4.kek.jp](https://bugzilla-geant4.kek.jp)

## The Geant4 Collaboration

Geant4 is developed and maintained by the Geant4 Collaboration, a worldwide team of physicists and software engineers.

For more information about the collaboration:
- Visit [https://cern.ch/geant4](https://cern.ch/geant4)

## License

All contributions to Geant4 are made under the Geant4 Software License.

See the [LICENSE](../../LICENSE) file for full license terms.

## Related Documentation

- [Build System](./build-system): Understanding the build process
- [Source Modules](./source-modules): Overview of code organization
- [Architecture](../architecture): Understanding Geant4's design
