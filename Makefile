# Makefile to generate directory structure and manage git operations

# Configuration
INCLUDE_EXTENSIONS := .md .ipynb .r .py .sh
EXCLUDE_FILES := update.py README.md
EXCLUDE_DIRS := .git src CV

# Generate directory structure and write to README.md
README.md:
	@echo "# Repository Directory Structure\n" > README.md
	@find . -type f \( $(foreach ext,$(INCLUDE_EXTENSIONS),-name "*$(ext)" -o) -false \) \
		$(foreach dir,$(EXCLUDE_DIRS),-not -path "./$(dir)/*") \
		$(foreach file,$(EXCLUDE_FILES),-not -name "$(file)") \
		-not -path "./README.md" | while read -r file; do \
			[ -n "$$file" ] && echo "Found: $$file" >&2; \
		done > /tmp/find_output
	@if [ -s /tmp/find_output ]; then \
		sort /tmp/find_output | while read -r file; do \
			dir=$$(dirname "$$file" | sed 's|^\./||'); \
			base=$$(basename "$$file"); \
			if [ "$$dir" = "." ]; then \
				echo "- [$$base]($$base)" >> README.md; \
			else \
				if ! grep -q "^- \*\*$$dir/\*\*" README.md; then \
					echo "- **$$dir/**" >> README.md; \
				fi; \
				echo "  - [$$base]($$base)" >> README.md; \
			fi; \
		done; \
	else \
		echo "No matching files found, keeping README.md with header only." >&2; \
	fi
	@rm -f /tmp/find_output
	@echo "Directory structure has been written to README.md"

# Commit changes to git
commit: README.md
	@git add .
	@git commit -m "Update directory structure: $$(date '+%Y-%m-%d %H:%M:%S')"
	@echo "Successfully committed changes"

# Push changes to remote repository
push: commit
	@git push
	@echo "Successfully pushed changes"

# Default target
all: README.md commit push

.PHONY: all commit push